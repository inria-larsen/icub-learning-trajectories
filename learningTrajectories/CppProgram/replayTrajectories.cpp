/* 
 * For information: oriane.dermy@inria.fr (11/07/16)
 * 
 * This program allow to replay the learned distribution that is send by matlab by the port /matlab/write.
 * To launch it :
 * 1. launch yarpserver
 * 2. launch gazebo (you can use the world "wolrdPROMPS 
 */

#include <cmath>
#include <string>
#include <algorithm>
#include <map>
#include <iostream>
#include <fstream>

#include <yarp/os/all.h>
#include <yarp/dev/all.h>
#include <yarp/sig/all.h>
#include <yarp/math/Math.h>

#include "cartesianClient.h"

#define DEG2RAD     (M_PI/180.0)
#define RAD2DEG     (180.0/M_PI)
#define MAX_TORSO_PITCH     0.0    // [deg]

using namespace std;
using namespace yarp::os;
using namespace yarp::sig;

class Test: public RFModule
{
protected:
    // cartesian
    CartesianClient client;
    Vector xd;
    Vector od;
    Vector fext;
    string part;
    
    //my data
    vector <vector <double> > fx, fy, fz; //Will be possible to learn other trajectory
    vector <double>  ftx, fty, ftz; //temporary vector
    int nbReplay;
	int nbTimeStep;
    BufferedPort<Bottle> port;
    int cpt=0;
    bool firstTime;
    
    bool initCartesian(const string &robot, const string &part, bool swap_x=false,bool swap_y=false)
    {                        
        client.init(robot,part,swap_x,swap_y); 
        port.open("/replay/read");
        int sys = system("yarp connect /matlab/write /replay/read");// we could verify it is equal to 0;
        sys = system("yarp connect /replay/read /matlab/write");// we could verify it is equal to 0;
        xd.resize(3);
        od.resize(4);
        fext.resize(3);
        od[0]=-0.039151; od[1]=0.601265; od[2]=-0.79809; od[3]=2.851949; //see serena-ivaldi/WoZ/app/basicfiles/action.conf for more information about this values
		fext[0] = 0.000001; fext[1] = 0.000001; fext[2] = 0.000001; //TODO put real data!!!
		firstTime=true;
        return true;
    }
    
    void closeCartesian()
    {    
      client.close();      
    }

public:
 
    Test(): RFModule()
    {
	
    }
    
    bool configure(ResourceFinder &rf)
    {    
		string name=rf.check("name",Value("test_feedback")).asString().c_str();
        string robot=rf.check("robot",Value("icubGazeboSim")).asString().c_str();
		string part=rf.check("part",Value("left_arm")).asString().c_str();	
        bool swap_x=rf.check("swap_x",Value("false")).asBool();
        bool swap_y=rf.check("swap_y",Value("false")).asBool();    
		
        if (!initCartesian(robot,part,1,1))
          return false;
        return true;
    }
    
    bool close()
    {
		cout << "Close the module." << endl;
		closeCartesian();
		cout << "Close the yarp port and connection." << endl;
		port.close();
		int sys = system("Yarp clean"); // we could verify it is equal to 0;
		return true;
    }

    double getPeriod()
    {	
        return 0.01;
    }

    bool updateModule()
    { 
       Vector buttons,rpy;
       Vector x,o;    
       //client.getPose(x,o); // get current pose as x,o

      //  cout << "waiting for input" << endl;
        Bottle *input = port.read();
        if (input!=NULL) {
            cout << "got " << input->toString().c_str() << endl;
            if(input ==0) cpt =0;
            else if(input->size() == 1){
				 if(input->get(0).asDouble() == -1.0)
				 {
					 cout << "Closing the programm" << endl;
					 return false;
				 }
				 else cout << "Error should receive 3 variables." << endl;
            }
            for (int i=0; i<3; i++) 
			{
                xd[i] = input->get(i).asDouble();
            }
            double dist;
            int nbIt=0;
            do
            {
				if(firstTime==true)
				{
					client.setTrajectoryTime(3.0);
					firstTime = false;
				}else
				{
					client.setTrajectoryTime(0.1);
				}
				client.goToPose(xd,od); // new target is xd,od
				client.getPose(x,o); // get current pose as x,o
				dist = (x[0] -xd[0])*(x[0]-xd[0]) + (x[1] -xd[1])*(x[1]-xd[1]) + (x[2] -xd[2])*(x[2]-xd[2]);
				nbIt++;
			}while(dist > 0.001 && nbIt< 100);
				yInfo("Current position = (%s)",x.toString(3,3).c_str());
				yInfo(" Target position = (%s)",xd.toString(3,3).c_str());
				yInfo("		 Max Forces = (%s)", fext.toString(3,3).c_str());
            Bottle& output = port.prepare();
            output.clear();
            //output.addString("OK: ");
            output.addDouble(fext[0]);
            output.addDouble(fext[1]);
            output.addDouble(fext[2]);
            //cout << "writing " << output.toString().c_str() << endl;
            port.write();
        }
        
        cpt++;
    
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

    Test test;    
    int r=test.runModule(rf);  
        
    return r;
}
