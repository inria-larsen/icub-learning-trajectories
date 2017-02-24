/* 
 * For information: oriane.dermy@inria.fr (11/07/16)
 * 
 * In this program, we can control the arm robot the simulated robot using the geomagic touch.
 * During all the movement the partner has to maintain the black button, to allow the robot to record data in the file text.txt. 
 * After that, if we press the white button, the robot will replay the movement. The geomagic touch arm robot will follow this movement.
 * 
 * Be careful, when you replay the trajectory, the geomagic movement is quite brutal, you should keep holding (in a compliant way) the pen of the geomagic to avoid big movement. 
 * 
 * Moreover the orientation of the arm is corrected to be natural.
 * 
 */

#include <cmath>
#include <string>
#include <algorithm>
#include <map>
#include <iostream>
#include <fstream>
#include <boost/lexical_cast.hpp>

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

using namespace std;
using namespace yarp::os;
using namespace yarp::dev;
using namespace yarp::sig;
//using namespace yarp::math;
using namespace hapticdevice;

class Test: public yarp::os::RFModule
{
protected:
    // geomagic parameters
    PolyDriver drvGeomagic;
    IHapticDevice *igeo;    
    bool fakehaptic;
    Vector maxFeedback;
    Vector feedback;
    double minForce;
    double maxForce;
    int verbosity;

    // cartesian parameters
    CartesianClient client;
    Vector xd;
    Vector od;
    string part;
    
    // my data
    
    vector <vector <double> > secTrajX, secTrajY, secTrajZ, secTime; 
    vector <double>  trajX, trajy, trajz, trajTime; 
    int nbReplay, totalTraj;
	int nbTimeStep;
    int flagReplay, flagLearn;
    ofstream record;
    BufferedPort<Bottle> port;
    double currentTime;
    string fileName;
        
    bool initFile()
    {
		fileName = string("record"+ boost::lexical_cast<string>(totalTraj) + ".txt");  
		record.open(fileName.c_str(), ios::out | ios::trunc); 
        if (!record)
        {  
			cerr << "Cannot creat the file " << totalTraj  << " for recording movements." << endl;
			return false;
		}else
		{
			cout << "Correctly created " << fileName << endl;
		}
		return true;
	}
	bool closeFile()
	{
			totalTraj++;
			record.close();
			cout << "in close file" << endl;
			return true;
	}





    bool initGeomagic(const string &name, const string &geomagic)
    {
        Property optGeo("(device hapticdeviceclient)");
        optGeo.put("remote",("/"+geomagic).c_str());
        optGeo.put("local",("/"+name+"/geomagic").c_str());
        if (!drvGeomagic.open(optGeo)) return false;
        drvGeomagic.view(igeo);              

        Matrix T=yarp::math::zeros(4,4);
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
		cout << "Closing geomagic..." << endl;	
		igeo->setTransformation(yarp::math::eye(4,4));
		igeo->stopFeedback();
        drvGeomagic.close();                
        cout << "Geomagic is closed." << endl;	
    }
    
    bool initCartesian(const string &robot, const string &part, bool swap_x=false,bool swap_y=false)
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
 
    Test(Network &yarp, int verbositylevel=2): RFModule()
    {

       port.open("/ori_record/read");  

	   yarp.connect("/wholeBodyDynamicsTree/left_arm/cartesianEndEffectorWrench:o","/ori_record/read", "tcp");//, false); 
    }
    
    bool configure(ResourceFinder &rf)
    {    
		string name=rf.check("name",Value("test_feedback")).asString().c_str();
        string robot=rf.check("robot",Value("icubGazeboSim")).asString().c_str();
		string part=rf.check("part",Value("left_arm")).asString().c_str();	
    	
        string geomagic=rf.check("geomagic",Value("geomagic")).asString().c_str();
        fakehaptic=rf.check("fakehaptic",Value("false")).asBool();
        minForce=fabs(rf.check("min-force-feedback",Value(0.01)).asDouble());
        maxForce=fabs(rf.check("max-force-feedback",Value(1.5)).asDouble());
        bool swap_x=rf.check("swap_x",Value("false")).asBool();
        bool swap_y=rf.check("swap_y",Value("false")).asBool();
                
        // my initiate values
        nbTimeStep = 0;
        flagReplay = 0;
        flagLearn = 0;
        nbReplay = 0;
        totalTraj = 0;

	
		// if we don't want to use the geomagic: not tested here
        if (fakehaptic)
        {
          igeo=new FakeGeomagicDriver();
          Matrix T=yarp::math::zeros(4,4);
          T(0,1)=1.0;
          T(1,2)=1.0;
          T(2,0)=1.0;
          T(3,3)=1.0;
          igeo->setTransformation(yarp::math::SE3inv(T));
          igeo->setCartesianForceMode();
        }
        else
        {  
			if (!initGeomagic(name,geomagic)) return false;
        }
        
        if (!initCartesian(robot,part,1,1)) return false;
        return true;
    }
    
    bool close()
    {
		cout << "Closing all modules" << endl;
		if (!fakehaptic)
		{
			closeGeomagic();
		}
		closeCartesian();
		cout << "All modules are closed" << endl;
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
        Vector buttons,pos,rpy;
        igeo->getButtons(buttons);
        igeo->getPosition(pos);
        igeo->getOrientation(rpy);
		
        bool b0=(buttons[0]!=0.0);
        bool b1=(buttons[1]!=0.0);

        xd=pos;
        od[0]=-0.039151; od[1]=0.601265; od[2]=-0.79809; od[3]=2.851949; //see serena-ivaldi/WoZ/app/basicfiles/action.conf for more information about this values
        
        Vector x,o;    
        client.getPose(x,o); // get current pose as x,o
        
        //Treat button to know in which phase we are
        if((buttons[0] != 0.0) && (flagLearn ==0))
        {
			flagLearn = 1;
			cout << "Replay one from " << secTrajX.size() << " trajectories" << endl;
			initFile();
		}
		else if((buttons[0] == 0.0) && (flagLearn == 1) && (trajX.size() > 2))
		{
			cout << "End of the recorded trajectory." << endl;
			secTrajX.push_back(trajX);
			secTrajY.push_back(trajy);
			secTrajZ.push_back(trajz);
			secTime.push_back(trajTime);
			cout << "We've learned " << secTrajX.size() << " trajectories." << endl;

			record << endl;
			closeFile();
			trajX.clear();
			trajy.clear();
			trajz.clear(); 
			trajTime.clear();
			flagLearn = 0;
		}
		if((buttons[1]!=0) && (flagReplay==0))
		{
			if (secTrajX.size() > 0)
			{
				cout << "Give the number of the trajectory you want to replay between 0 and " << secTrajX.size() -1 << endl;
				cin >> nbReplay;
				if( (nbReplay < 0) || nbReplay > secTrajX.size() - 1)
				{
					cout << "I don't understand your answer, so I replay the first one." << endl;
					nbReplay = 0;
				}
			
				flagReplay = 1;
				nbTimeStep = 0;
			 }
			 else cout << "I cannot replay since I haven't learned any movements." << endl;	 
		}
		
		//if we learn
		if(flagLearn == 1)
		{
			/*If you don't want to record forces delet these 6 lines*/
			Bottle *input = port.read();
			if (input!=NULL) 
			{
				if(verbosity == 1) cout << "Got: " << input->toString().c_str() << endl;
			}
			else cout << "Got nothing" << endl;
		
			currentTime = Time::now();
			trajTime.push_back(currentTime);
			trajX.push_back(x[0]);
			trajy.push_back(x[1]);
			trajz.push_back(x[2]);
			cout << "Learn = " << trajX.size() << endl;
			record << currentTime << " " << x[0] << " " << x[1] << " " << x[2]  << " " << input->get(0).asDouble() << " " << input->get(1).asDouble() << " " << input->get(2).asDouble() << endl;
			cout << "record ..." << endl;
			client.goToPose(xd,od); // new target is xd,od
			//feedback=yarp::math::operator*(yarp::math::operator-(x,xd), 45.0);	
			 feedback=yarp::math::operator*(x,0.0);	
			verifyFeedback();
			igeo->setFeedback(feedback);
		}
		
		else if(flagReplay == 1) // if the button replay is pushed
		{
			cout << "Replay. Time step = " << nbTimeStep << endl;

			xd[0] = secTrajX[nbReplay][nbTimeStep];
			xd[1] = secTrajY[nbReplay][nbTimeStep];
			xd[2] = secTrajZ[nbReplay][nbTimeStep];	
			nbTimeStep++;
			
			// the sync don't work with icubSim
			client.goToPose(xd,od); // new target is xd,od
			//client.waitMotionDone(0.04);
			client.getPose(x,o);
			igeo->getPosition(pos);
			double distance= ((x[0] - pos[0])*(x[0] - pos[0]) + (x[1] - pos[1])*(x[1] - pos[1]) + (x[2] - pos[2])*(x[2] - pos[2]));
			while(distance > 0.001)
			{
				cout << "Distance =" << distance << endl;
				yInfo("Sim position    = (%s)",x.toString(3,3).c_str());
				yInfo("Haptic position = (%s)", pos.toString(3,3).c_str());
				feedback=yarp::math::operator*(yarp::math::operator-(x,pos), 45.0); //(x-xd)*45.0; 

				
				yInfo("Feedback = (%s)", feedback.toString(3,3).c_str());
				verifyFeedback(); 
				yInfo("Feedback = (%s)", feedback.toString(3,3).c_str());

				igeo->setFeedback(feedback);			
				igeo->getPosition(pos);
				client.getPose(x,o);
				distance= ((x[0] - pos[0])*(x[0] - pos[0]) + (x[1] - pos[1])*(x[1] - pos[1]) + (x[2] - pos[2])*(x[2] - pos[2]));
			}
			if (nbTimeStep >= secTrajX[0].size())
			{
				cout << "End of the trajectory." << endl;
				nbTimeStep = 0;
				flagReplay = 0;
			}
		}
		else
		{
			client.goToPose(xd,od); // new target is xd,od
			feedback=yarp::math::operator*(x,0.0);	
			verifyFeedback();
			igeo->setFeedback(feedback);
		/** if we want to force the icub and the geomagic to take the origin position do:.
		 * //we force the icub to take the origin position 
		*		//client.goToPose(xd,od); // new target is xd,od
		*		//feedback=yarp::math::operator*(x,0.0);	
		*		//verifyFeedback();
		*	    //igeo->setFeedback(feedback);
		*		xd[0] = -0.04939;
		*		xd[1] = -0.06839;
		*		xd[2] = -0.05287;
		*		double distance;
		*		do
		*		{
		*			client.goToPose(xd,od); // new target is xd,od
		*			client.getPose(x,o);
		*			distance= ((x[0] - xd[0])*(x[0] - xd[0]) + (x[1] - xd[1])*(x[1] - xd[1]) + (x[2] - xd[2])*(x[2] - xd[2]));
		*		}while(distance > 0.001);
		*		
		*		//we force the geom to take the origin position
		*		igeo->getPosition(pos);
		*		distance= ((xd[0] - pos[0])*(xd[0] - pos[0]) + (xd[1] - pos[1])*(xd[1] - pos[1]) + (xd[2] - pos[2])*(xd[2] - pos[2]));
		*		//cout << "distance init= " << distance << " x= " << x.toString(3,3).c_str() << "pos=" <<  pos.toString(3,3).c_str()<< endl;
		*		while(distance > 0.00005)
		*		{
		*			cout << "positionning the geomagic... =" << distance << endl;
		*		    //yInfo("Sim position    = (%s)",x.toString(3,3).c_str());
		*			//yInfo("Haptic position = (%s)", pos.toString(3,3).c_str()); 
		*			feedback=yarp::math::operator*(yarp::math::operator-(xd,pos), 50.0); //(x-xd)*45.0; 
		*
		*			//yInfo("Feedback = (%s)", feedback.toString(3,3).c_str());
		*			verifyFeedback(); 
		*			//yInfo("Feedback = (%s)", feedback.toString(3,3).c_str());
	    *
		*			igeo->setFeedback(feedback);			
		*			igeo->getPosition(pos);
		*			client.getPose(x,o);
		*			distance= ((xd[0] - pos[0])*(xd[0] - pos[0]) + (xd[1] - pos[1])*(xd[1] - pos[1]) + (xd[2] - pos[2])*(xd[2] - pos[2]));		
		*		}
		*		cout << "Geomagic positionned" << endl;
		***/
		
				
		}
		     
        return true;
    }    
};



int main(int argc,char *argv[])
{           
    Network yarp;

    if (!yarp.checkNetwork())
    {
        yError("YARP server not found!");
        return 1;
    }

    ResourceFinder rf;
    rf.configure(argc,argv);

    Test test(yarp);    
    int r=test.runModule(rf);  
    
    return r;
}
