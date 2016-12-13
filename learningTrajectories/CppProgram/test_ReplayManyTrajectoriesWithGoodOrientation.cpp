/* 
 * For information: oriane.dermy@inria.fr (11/07/16)
 * 
 * This program is based on the main.cpp program which is an example of the utilization of the arm robot "geomagic Touch".
 * Here, the simulated robot follows the trajectory given by the partner on the geomagic touch arm robot.
 * During all the movement the partner has to maintain the black button, to allow the robot to learn it. After that, if we press the white button, the robot will replay the movement. The geomagic touch arm robot will try to follow this movement.
 * 
 * Moreother the orientation of the arm is corrected to be natural.
 * LIMITATION: the arm robot follow the movement in a hard way (not smooth at all).
 * 
 * 
 */

#include <cmath>
#include <string>
#include <algorithm>
#include <map>

#include <yarp/os/all.h>
#include <yarp/dev/all.h>
#include <yarp/sig/all.h>
#include <yarp/math/Math.h>

#include <hapticdevice/IHapticDevice.h>

#include "fakeGeomagicDriver.h"
#include "cartesianClient.h"

#define DEG2RAD     (M_PI/180.0)
#define RAD2DEG     (180.0/M_PI)
#define MAX_TORSO_PITCH     0.0    // [deg]

//using namespace std;
//using namespace yarp::os;
//using namespace yarp::dev;
//using namespace yarp::sig;
//using namespace yarp::math;
using namespace hapticdevice;

class Test: public yarp::os::RFModule
{
protected:
    // geomagic  
    yarp::dev::PolyDriver drvGeomagic;
    IHapticDevice *igeo;    
    bool fakehaptic;
    yarp::sig::Vector maxFeedback;
    yarp::sig::Vector feedback;
    double minForce;
    double maxForce;
    int verbosity;

    // cartesian
    CartesianClient client;
    yarp::sig::Vector xd;
    yarp::sig::Vector od;
    
    std::string part;
    
    //my data
    std::vector <std::vector <double> > fx, fy, fz; //Will be possible to learn other trajectory
    std::vector <double>  ftx, fty, ftz; //temporary vector
    int nbReplay;
	int nbTimeStep;
    int flagReplay, flagLearn;
    
    bool initGeomagic(const std::string &name, const std::string &geomagic)
    {
        yarp::os::Property optGeo("(device hapticdeviceclient)");
        optGeo.put("remote",("/"+geomagic).c_str());
        optGeo.put("local",("/"+name+"/geomagic").c_str());
        if (!drvGeomagic.open(optGeo))
            return false;
        drvGeomagic.view(igeo);              

        yarp::sig::Matrix T=yarp::math::zeros(4,4);
        T(0,1)=1.0;
        T(1,2)=1.0;
        T(2,0)=1.0;
        T(3,3)=1.0;
        igeo->setTransformation(yarp::math::SE3inv(T));
        igeo->setCartesianForceMode();                
        feedback.resize(3,0.0);
        return true;
    }
    
    void closeGeomagic()
    {
		if (verbosity >=1 ) std::cout << "closing geomagic " << std::endl;
		
		if (verbosity >= 2) std::cout << "..... set transf" << std::endl;
		igeo->setTransformation(yarp::math::eye(4,4));
		

		if (verbosity >= 2) std::cout << "..... stop feedback" << std::endl;
		igeo->stopFeedback();
		
		
        if (verbosity >= 2) std::cout << "..... stop feedback" << std::endl;
        drvGeomagic.close();                
        
    }
    
    bool initCartesian(const std::string &robot, const std::string &part, bool swap_x=false,bool swap_y=false)
    {                        
        client.init(robot,part,swap_x,swap_y);
       
        xd.resize(3);
        od.resize(4);
        return true;
    }
    
    void closeCartesian()
    {    
      client.close();      
    }

public:
 
    Test(int verbositylevel=2): RFModule()
    {
		verbosity = verbositylevel;
		if (verbosity == 2) std::cout << "Verbosity set to "<<verbosity<<std::endl;
    }
    
    bool configure(yarp::os::ResourceFinder &rf)
    {    
        std::string name=rf.check("name",yarp::os::Value("test_feedback")).asString().c_str();
        std::string robot=rf.check("robot",yarp::os::Value("icubGazeboSim")).asString().c_str();
		std::string part=rf.check("part",yarp::os::Value("left_arm")).asString().c_str();	
    	
        std::string geomagic=rf.check("geomagic",yarp::os::Value("geomagic")).asString().c_str();
        fakehaptic=rf.check("fakehaptic",yarp::os::Value("false")).asBool();
        minForce=fabs(rf.check("min-force-feedback",yarp::os::Value(0.01)).asDouble());
        maxForce=fabs(rf.check("max-force-feedback",yarp::os::Value(1.5)).asDouble());
        bool swap_x=rf.check("swap_x",yarp::os::Value("false")).asBool();
        bool swap_y=rf.check("swap_y",yarp::os::Value("false")).asBool();
                
        // my initiate values
        nbTimeStep = 0;
        flagReplay = 0;
        flagLearn = 0;
        nbReplay = 0;
        
        if (fakehaptic)
        {
          igeo=new FakeGeomagicDriver();
          yarp::sig::Matrix T=yarp::math::zeros(4,4);
          T(0,1)=1.0;
          T(1,2)=1.0;
          T(2,0)=1.0;
          T(3,3)=1.0;
          igeo->setTransformation(yarp::math::SE3inv(T));
          igeo->setCartesianForceMode();
        }
        else
        {  
        if (!initGeomagic(name,geomagic))
          return false;
        }
        if (!initCartesian(robot,part,1,1))
          return false;
        return true;
    }
    
    bool close()
    {
		std::cout << "close the modul" << std::endl;
      if (!fakehaptic)
      {
        closeGeomagic();
      }
      closeCartesian();
      return true;
    }

    double getPeriod()
    {	
        return 0.01;
    }


	void verifyFeedback()                     
	{
        for (size_t i=0;i<3;i++)
        {
          if (feedback[i]>maxForce) feedback[i]=maxForce;
          if (feedback[i]<-maxForce) feedback[i]=-maxForce;
          if ((feedback[i]>=0) && (feedback[i]<minForce)) feedback[i]=0.0; // remove small values
          if ((feedback[i]<0) && (feedback[i]>-minForce)) feedback[i]=0.0; // remove small values
        }
	}

    bool updateModule()
    { 
        yarp::sig::Vector buttons,pos,rpy;
        igeo->getButtons(buttons);
        igeo->getPosition(pos);
        igeo->getOrientation(rpy);
		
        bool b0=(buttons[0]!=0.0);
        bool b1=(buttons[1]!=0.0);

        xd=pos;
        od[0]=-0.039151; od[1]=0.601265; od[2]=-0.79809; od[3]=2.851949; //see serena-ivaldi/WoZ/app/basicfiles/action.conf for more information about this values
        
        yarp::sig::Vector x,o;    
        client.getPose(x,o); // get current pose as x,o
        
        //Treat button to know in which phase we are
        if((buttons[0] != 0.0) && (flagLearn ==0))
        {
			flagLearn = 1;
			std::cout << "Replay one from " << fx.size() << " trajectories" << std::endl;
		}
		else if((buttons[0] == 0.0) && (flagLearn == 1) && (ftx.size() > 2))
		{
			std::cout << "end of the learned trajectory" << std::endl;
			fx.push_back(ftx);
			fy.push_back(fty);
			fz.push_back(ftz);
			std::cout << "We've learned " << fx.size() << " trajectories" << std::endl;

			ftx.clear();
			fty.clear();
			ftz.clear(); 
			flagLearn = 0;
		}
		if((buttons[1]!=0) && (flagReplay==0))
		{
			if (fx.size() > 0)
			{
				std::cout << "give the number of replay between 0 and " << fx.size() -1 << std::endl;
				std::cin >> nbReplay;
				if( (nbReplay < 0) || nbReplay > fx.size() - 1)
				{
					std::cout << "error of value, we replay the first one" << std::endl;
					nbReplay = 0;
				}
			
				flagReplay = 1;
				nbTimeStep = 0;
			 }
			 else std::cout << "We cannot replay since we don't learn any movements." << std::endl;
			 
		}
		
		//if we learn
		if(flagLearn == 1)
		{
			ftx.push_back(x[0]);
			fty.push_back(x[1]);
			ftz.push_back(x[2]);
			std::cout << "learn = " << ftx.size() << std::endl;
			client.goToPose(xd,od); // new target is xd,od
			feedback=yarp::math::operator*(yarp::math::operator-(x,xd), 45.0);
			
			 
			verifyFeedback();
			igeo->setFeedback(feedback);
		}
		
		else if(flagReplay == 1) // if the button replay is pushed
		{
			
			std::cout << "replay. Time step = " << nbTimeStep << std::endl;

			xd[0] = fx[nbReplay][nbTimeStep];
			xd[1] = fy[nbReplay][nbTimeStep];
			xd[2] = fz[nbReplay][nbTimeStep];	
			nbTimeStep++;
			
			client.goToPose(xd,od); // new target is xd,od
			client.getPose(x,o);
			igeo->getPosition(pos);
			int wait;
			double distance= ((x[0] - pos[0])*(x[0] - pos[0]) + (x[1] - pos[1])*(x[1] - pos[1]) + (x[2] - pos[2])*(x[2] - pos[2]));
			while(distance > 0.001)
			{
				
				std::cout << "distance =" << distance << std::endl;
			    yInfo("Sim position    = (%s)",x.toString(3,3).c_str());
				yInfo("Haptic position = (%s)", pos.toString(3,3).c_str());
				feedback=yarp::math::operator*(yarp::math::operator-(x,pos), 45.0); //(x-xd)*45.0; 
				
		//		std::cout<< "pos -x "<< 3.0*pos  << std::endl;
				
				//feedback=yarp::math::operator*(yarp::math::operator-(x,xd), 20.0); //(x-pos)*45.0; // the feedback force will be proportional to the error in position    
				/*std::cout << "val : " << abs(abs(x[0]) - abs(pos[0])) << std::endl;
				feedback[0] = abs(abs(x[0]) - abs(pos[0]));
				feedback[1] = abs(abs(x[1]) - abs(pos[1]));
				feedback[2] = abs(abs(x[2]) - abs(pos[2]));
				if(x[0] < pos[0]) feedback[0] -=  feedback[0];
				if(x[1] < pos[1]) feedback[1] -=  feedback[1];
				if(x[2] < pos[2]) feedback[2] -=  feedback[2];*/
				
				yInfo("Feedback              = (%s)", feedback.toString(3,3).c_str());
				verifyFeedback(); 
				yInfo("after verify Feedback = (%s)", feedback.toString(3,3).c_str());

			//	yInfo("        feedback = (%s)",feedback.toString(3,3).c_str());
				igeo->setFeedback(feedback);			
				igeo->getPosition(pos);
				client.getPose(x,o);
				distance= ((x[0] - pos[0])*(x[0] - pos[0]) + (x[1] - pos[1])*(x[1] - pos[1]) + (x[2] - pos[2])*(x[2] - pos[2]));
			
				
			}
	
			
			if (nbTimeStep >= fx[0].size())
			{
				std::cout << "end of the trajectory" << std::endl;
				nbTimeStep = 0;
				flagReplay = 0;
			}
	
		}
		else
		{
			client.goToPose(xd,od); // new target is xd,od
			feedback=yarp::math::operator*(yarp::math::operator-(x,xd), 45.0); //(x-xd)*45.0; 
			
			
			verifyFeedback();
			igeo->setFeedback(feedback);
		 //client.getPose(x,o); // get current pose as x,o

			//yInfo("current position = (%s)",x.toString(3,3).c_str());
			// yInfo(" target position = (%s)",xd.toString(3,3).c_str());
			// yInfo("        feedback = (%s)",feedback.toString(3,3).c_str());
			 //yInfo("    Position arm = (%s)", pos.toString(3,3).c_str());
		}
		


        
        return true;
    }    
};



int main(int argc,char *argv[])
{           
    yarp::os::Network yarp;

    if (!yarp.checkNetwork())
    {
        yError("YARP server not found!");
        return 1;
    }

    yarp::os::ResourceFinder rf;
    rf.configure(argc,argv);

    Test test;    
    int r=test.runModule(rf);  
    
    return r;
}
