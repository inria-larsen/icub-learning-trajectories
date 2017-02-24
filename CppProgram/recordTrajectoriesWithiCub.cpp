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
 * Thanks to this program, we can record the trajectory and the forces of the left arm of the icub robot. 
 * These data are recorded in the file record.txt, when we press the button "s" followed by a chariot return.
 * To stop the recording, press again s.
 * If you want to visualize the movement, press the button "r" (replay).
 * 
 * 
 * Moreother the orientation of the arm is corrected to be natural.
 * 
 */
 
 
 /// WARNING: THIS FUNCTION IS IN DEVELOPMENT, DO NOT USE IT WITHOUT VERIFY
 

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
using namespace yarp::dev;
using namespace yarp::sig;


//a function I found on the net to be not blocked during the lecture of data.
int kbhit(void)
{
  struct timeval tv;
  fd_set read_fd;

  /* Do not wait at all, not even a microsecond */
  tv.tv_sec=0;
  tv.tv_usec=0;

  /* Must be done first to initialize read_fd */
  FD_ZERO(&read_fd);

  /* Makes select() ask if input is ready:
   * 0 is the file descriptor for stdin    */
  FD_SET(0,&read_fd);

  /* The first parameter is the number of the
   * largest file descriptor to check + 1. */
  if(select(1, &read_fd,NULL, /*No writes*/NULL, /*No exceptions*/&tv) == -1)
    return 0;  /* An error occured */

  /*  read_fd now holds a bit map of files that are
   * readable. We test the entry for the standard
   * input (file 0). */
  
if(FD_ISSET(0,&read_fd))
    /* Character pending on stdin */
    return 1;

  /* no characters were pending */
  return 0;
}

class Test: public RFModule
{
protected:

    // cartesian
    CartesianClient client;
    Vector xd, xinit;
    Vector od;
    string part;
    
    //my data
    vector <vector <double> > fx, fy, fz; //Will be possible to learn other trajectory
    vector <double>  ftx, fty, ftz; //temporary vector
    int nbReplay;
	int nbTimeStep;
    int flagReplay, flagLearn;
    ofstream record;
    int verbositylevel;
	char c;
	BufferedPort<Bottle> port;
	
    bool initCartesian(const string &robot, const string &part, bool swap_x=false,bool swap_y=false)
    {      		
        client.init(robot,part,swap_x,swap_y);
        xd.resize(3);
        xinit.resize(3);
        od.resize(4);
        client.getPose(xinit,od);
        od[0]=-0.039151; od[1]=0.601265; od[2]=-0.79809; od[3]=2.851949; //see serena-ivaldi/WoZ/app/basicfiles/action.conf for more information about this values

        return true;
    }
    
    void closeCartesian()
    {    
      client.close();      
    }

public:
 
    Test(Network &yarp, int verbose=1): RFModule()
    {
        verbositylevel = verbose; 
        port.open("/ori_record/read");  
        //TODO to verify it is this port that give information about external forces.
	    yarp.connect("/wholeBodyDynamics/left_arm/ext_ft_sens:o","/ori_record/read", "tcp");//, false);      
    }
    
    bool configure(ResourceFinder &rf)
    {    
        string robot=rf.check("robot",Value("icub")).asString().c_str();
		string part=rf.check("part",Value("left_arm")).asString().c_str();	

        bool swap_x=rf.check("swap_x",Value("false")).asBool();
        bool swap_y=rf.check("swap_y",Value("false")).asBool();
                 
        // my initiate values
        nbTimeStep = 0;
        flagReplay = 0;
        flagLearn = 0;
        nbReplay = 0;
		record.open("record.txt", ios::out | ios::trunc); 
        if (!record)  cerr << "Cannot creat the file for recording movements." << endl;

        if (!initCartesian(robot,part,1,1)) return false;
        return true;
    }
    
    bool close()
    {
		if(verbositylevel == 1) cout << "Closing all modules ... ";
		record.close();
		port.close();
		if(verbositylevel == 1) cout << " ... All modules are closed." << endl;
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

   //     client.getPose(x,o); // get current pose as x,o
        
		if( kbhit() != 0 ) //Keys interaction
		{        
			c = getchar();
			if(verbositylevel == 1) cout << "button pressed: " <<  c << endl;
            
            if((c=='s') && (flagLearn == 0)) //begin to record
            {
				flagLearn = 1;		
				if(verbositylevel == 1)  cout << "Record " << fx.size()+1 << " trajectories" << endl;
			}else if((c == 's') && (flagLearn == 1) && (ftx.size() > 2)) //stop recording
			{
				cout << "End of the recorded trajectory." << endl;
				client.setTrajectoryTime(3.0);
				client.goToPose(xinit,od); // new target is xd,od
				client.setTrajectoryTime(0.5);
				if(verbositylevel == 1) cout << "Rythme slow to replace the arm." << endl;
				fx.push_back(ftx);
				fy.push_back(fty);
				fz.push_back(ftz);
				if(verbositylevel == 1) cout << "We've learned " << fx.size() << " trajectories." << endl;
	
				record << endl;
				ftx.clear();
				fty.clear();
				ftz.clear(); 
				flagLearn = 0;
			}else if((c== 'r')  && (flagReplay==0))
			{
				if (fx.size() > 0)
				{
					cout << "Give the number of the trajectory you want to replay between 0 and " << fx.size() -1 << endl;
					cin >> nbReplay;
					if( (nbReplay < 0) || nbReplay > fx.size() - 1)
					{
						cout << "I don't understand your answer, so I replay the first one." << endl;
						nbReplay = 0;
					}
			
					flagReplay = 1;
					nbTimeStep = 0;
				}
				else cout << "I cannot replay since I haven't learned any movements." << endl;	 
			}    
        }
	
		//if you are recording
		if(flagLearn == 1)
		{
			client.getPose(x,o); // get current pose as x,o
			
			Bottle *input = port.read();
			if (input!=NULL) 
			{
				if(verbositylevel == 1) cout << "Got: " << input->toString().c_str() << endl;
			}
			else cout << "Got nothing" << endl;
		
			ftx.push_back(x[0]);
			fty.push_back(x[1]);
			ftz.push_back(x[2]);
			cout << "It records: " << ftx.size() << endl;
			record << x[0] << " " << x[1] << " " << x[2] << " " << input->get(0).asDouble() << " " << input->get(1).asDouble() << " " << input->get(2).asDouble() << endl;		
			xd = x;
			//client.goToPose(xd,od); // new target is xd,od
		}
		
		else if(flagReplay == 1) // if you are replaying
		{
			
			cout << "Replay. Time step = " << nbTimeStep << endl;

			xd[0] = fx[nbReplay][nbTimeStep];
			xd[1] = fy[nbReplay][nbTimeStep];
			xd[2] = fz[nbReplay][nbTimeStep];	
			nbTimeStep++;
			
			client.goToPose(xd,od); 
			client.getPose(x,o);
		
			yInfo("Computed position = (%s)",x.toString(3,3).c_str());
			yInfo("Desired position  = (%s)",xd.toString(3,3).c_str());

			if (nbTimeStep >= fx[0].size())
			{
				cout << "End of the trajectory." << endl;
				if(verbositylevel == 1) cout << "Rythme slow to replace the arm." << endl;

				client.setTrajectoryTime(3.0);
				//client.goToPoseSync(xinit,od); // new target is xd,od
				client.goToPose(xinit,od);   // send request and wait for reply
			//	client.waitMotionDone(0.04);
				client.setTrajectoryTime(0.5);
				nbTimeStep = 0;
				flagReplay = 0;
			}
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
