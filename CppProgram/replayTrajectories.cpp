
/* 
 * Copyright (C) <2016> RobotCub Consortium, European Commission FP6 Project IST-004370
 * Author Oriane Dermy
 * email:   oriane.dermy@inria.fr
 * website: www.robotcub.org
 * Permission is granted to copy, distribute, and/or modify this program
 * under the terms of the GNU General Public License, version 2 or any
 * later version published by the Free Software Foundation.
 *
 * A copy of the license can be found at
 * http://www.robotcub.org/icub/license/gpl.txt
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details
 * 
 * This program allow to replay the learned distribution that is send by matlab by the port /matlab/write.
 * The algorithm is not finished yet, I have to add the compliance of the arm.
 * 
 * To launch it :
 * 1. launch yarpserver
 * 2. launch iKin for the left arm
 * 2. launch wholebodydynamics (to have information about forces)
 * 4. launch the gravity compensator (to compensate gravity)
 * 2. launch the program on matlab
 * 3. launch this program.
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


///**********************************/
//impedance_velocity_mode off
//impedance_stiffness 0.5 0.5 0.5 0.2 0.1
//impedance_damping 60.0 60.0 60.0 20.0 0.0


 //impStiff[i]
   //impDamp[i]

    //IPositionControl  *posArm;
    //Vector leftArmReachOffs;
    //Vector leftArmGraspOffs;
    //Vector leftArmGraspSigma;
    //Vector leftArmHandOrien;
    //Vector leftArmJointsStiffness;
    //Vector leftArmJointsDamping;


        //IControlMode2    *imode=modeArm;
        //IPositionControl *ipos=posArm;
        //yInfo("*** Homing %s",type.c_str());
        //for (size_t j=0; j<homeVels.length(); j++)
            //imode->setControlMode(j,VOCAB_CM_POSITION);

        //for (size_t j=0; j<homeVels.length(); j++)
        //{
            //ipos->setRefSpeed(j,homeVels[j]);
            //ipos->positionMove(j,homePoss[j]);
        //}
        
        
                //yInfo("*** Stopping %s joints",type.c_str());
        //for (size_t j=0; j<homeVels.length(); j++)
            //imode->setControlMode(j,VOCAB_CM_POSITION);

        //for (size_t j=0; j<homeVels.length(); j++)
        //{
            //double fb;
            //ienc->getEncoder(j,&fb);
            //ipos->positionMove(j,fb);
        //}
    //}
    
    //// ouvrir ou fermer la main
            //yInfo("*** %s %s",actionStr.c_str(),type.c_str());
        //for (size_t j=0; j<handVels.length(); j++)
            //imode->setControlMode(homeVels.length()+j,VOCAB_CM_POSITION);

        //for (size_t j=0; j<handVels.length(); j++)
        //{
            //int k=homeVels.length()+j;
            //ipos->setRefSpeed(k,handVels[j]);
            //ipos->positionMove(k,(*poss)[j]);
        //}
        
                //// home part
        //Bottle &bHome=rf.findGroup("home_arm");
        //bHome.setMonitor(rf.getMonitor());
        //homePoss.resize(7,0.0); homeVels.resize(7,0.0);
        //if (!getHomeOptions(bHome, homePoss, homeVels)) { yError ("Error in parameters section 'home_arm'"); return false; }
    
    
    
           //if (useLeftArm)
        //{
            //drvCartLeftArm=new PolyDriver;
            //if (!drvCartLeftArm->open(optCartLeftArm))
            //{
                //close();
                //return false;
            //}

            //if (leftArmImpVelMode)
            //{
                //IInteractionMode  *imode;
                //IImpedanceControl *iimp;

                //drvLeftArm->view(imode);
                //drvLeftArm->view(iimp);

                //int len=leftArmJointsStiffness.length()<leftArmJointsDamping.length()?
                        //leftArmJointsStiffness.length():leftArmJointsDamping.length();

                //for (int j=0; j<len; j++)
                //{
                    //iimp->setImpedance(j,leftArmJointsStiffness[j],leftArmJointsDamping[j]);
                    //imode->setInteractionMode(j,VOCAB_IM_COMPLIANT);
                //}
            //}
        //}

    
    
///**********************************/


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
    BufferedPort<Bottle> port, portForces;
    int cpt=0;
    double compliance;
    bool firstTime;
    int verbositylevel;
    
    bool initCartesian(const string &robot, const string &part, bool swap_x=false,bool swap_y=false)
    {                        
        client.init(robot,part,swap_x,swap_y); 
        client.setPosePriority("position");
        xd.resize(3);
        
        od.resize(4);
        fext.resize(3);
        od[0]=-0.039151; od[1]=0.601265; od[2]=-0.79809; od[3]=2.851949; //see serena-ivaldi/WoZ/app/basicfiles/action.conf for more information about this values
		//fext[0] = 0.00048941; fext[1] =  0.00048941; fext[2] =  0.00048941; //TODO put real data!!!
		firstTime=true;
		compliance = 1;
        return true;
    }
    
    void closeCartesian()
    {    
      client.close();      
    }

public:
 
    Test(Network &yarp, int verbose = 1): RFModule()
    {
	    port.open("/replay/read");  //to comunicate with the matlab program.
	    portForces.open("/replay/readForces");   // to read forces from the wholebodyDynamics program.
	    yarp.connect("/matlab/write","/replay/read", "tcp");//, false);
	    yarp.connect("/replay/read", "/matlab/write","tcp");
	    yarp.connect("/wholeBodyDynamics/left_arm/ext_ft_sens:o","/replay/readForces", "tcp");

	    verbositylevel = verbose; 
    }
    
    bool configure(ResourceFinder &rf)
    {    
		string name=rf.check("name",Value("test_feedback")).asString().c_str();
        string robot=rf.check("robot",Value("icubNancy01")).asString().c_str();
		string part=rf.check("part",Value("left_arm")).asString().c_str();	
        bool swap_x=rf.check("swap_x",Value("false")).asBool();
        bool swap_y=rf.check("swap_y",Value("false")).asBool();    
		
        if (!initCartesian(robot,part,1,1))
          return false;
        return true;
    }
    
    bool close()
    {
		if(verbositylevel == 1) cout << "Close the module." << endl;
		closeCartesian();
		if(verbositylevel == 1) cout << "Close the yarp port and connection." << endl;
		port.close();
		portForces.close();
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

		if(verbositylevel == 1) cout << "Before reading" << endl;
        Bottle *input = port.read();
        if(verbositylevel == 1) cout << "Before reading forces." << endl;
        //Bottle *inputForces = portForces.read();
        if (input!=NULL) 
        {
            if(verbositylevel == 1) cout << "Got: " << input->toString().c_str() << endl;
            if(input->size() == 1)
            {
				 if(input->get(0).asDouble() == -1.0)
				 {
					 if(verbositylevel == 1) cout << "Closing the programm." << endl;
					 return false;
				 }
				 else if(input->get(0).asDouble() == 0.0)
				 {
					 if(verbositylevel == 1) cout << "New movement." << endl;
					 firstTime = true;
					 cpt=0;
					 Bottle& output = port.prepare();
					 output.clear();
					 output.addDouble(0);
					 port.write();
					 return true;
				 }
            }
            for (int i=0; i<3; i++) 
            {
                xd[i] = input->get(i).asDouble();
                //if(inputForces !=NULL)
                //{
					//fext[i] = inputForces->get(i).asDouble();
				//}else
				//{
					fext[i] = 0.00048941+ cpt*0.00001;
				//}
            }
            compliance = input->get(3).asDouble();
            double dist;
            int nbIt=0;
            do
            {
				if(firstTime==true)
				{
					client.setTrajectoryTime(3.0);
					if(verbositylevel == 1) cout << "Rythme slow to begin the movement" << endl;
					//doesn't work with icubGazebo
					//client.goToPoseSync(xd,od);   // send request and wait for reply
					client.goToPose(xd,od);
					//client.waitMotionDone(0.04);  // wait until the motion is done and ping at each 0.04 seconds
					firstTime = false;
				}else
				{
					client.setTrajectoryTime(0.5);
					client.goToPose(xd,od); // new target is xd,od
					//client.waitMotionDone(0.004);
				}
				
				client.getPose(x,o); // get current pose as x,o
				dist = (x[0] -xd[0])*(x[0]-xd[0]) + (x[1] -xd[1])*(x[1]-xd[1]) + (x[2] -xd[2])*(x[2]-xd[2]);
				nbIt++;
			}while(dist > 0.001 && nbIt< 100);
			
			//yInfo("Current position = (%s)",x.toString(3,3).c_str());
			yInfo(" Target position = (%s)",xd.toString(3,3).c_str());
			yInfo(" Compliance      = (%f)",compliance);
			
			//if(verbositylevel == 1)
			//{
				//bool value;
				//client.getReferenceMode(&value);
				//if(value == true) yInfo(" reference mode = true");
				//else yInfo(" reference mode = false");
			//}
            Bottle& output = port.prepare();
            output.clear();
            output.addDouble(fext[0]);
            output.addDouble(fext[1]);
            output.addDouble(fext[2]);
            if(verbositylevel == 1) cout << "It is writing: " << output.toString().c_str() << endl;
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

    Test test(yarp);    
    int r=test.runModule(rf);  
        
    return r;
}
