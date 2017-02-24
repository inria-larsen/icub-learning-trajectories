# icub-learning-trajectories
[under dev] One C++ program allow to record trajectories (forces + end effector's cartesian positions). Then, a matlab program computes the distribution other three kind of trajectories and allows to recognize and to continue a movement when it is initate by the icub's partner.

If you have some question, remarks etc, send an email to oriane.dermy@inria.fr.

***WARNING: this program is under developpment, do not launch it now without verify it.***

## PRE-INSTALLATION:
You need to have installed the geomagic touch (and of course gazebo, yarp, and so on).

## INSTALLATION:
`cd CppProgram`   
`mkdir build`   
`cd build`   
`ccmake ../`   
`make`   

REMARK:
In the matlab codes, the figure that represent the learned distribution and other interesting plot are commented. If you want to see this plot, you have to decomment them in the matlab code.

## V1 information
In this version, in the .zip file:   

A. you can launch a program that learn trajectories from the geomagic touch (see record trajectories.cpp). It will record a "record.txt" file. For example, you can use the world "worldPROMPS.sdf" to have some goal to achieve with the robot left arm. Open the code of this function to have more information about how to launch it.
It requires to have installed the geomagic touch driver.

B. You can learn the distribution of your trajectories. For example:  
1. Using the previous program, create three files of trajectories (for each file, you can do several samples) and replace the files in MatlabProgram.Data by them.   
2. Launch a yarpserver.   
3. Launch proMPs_withPhasis.m on Matlab. It will compute the distribution of your trajectories. This programm will wait with the message "Please connect to a bottle sink (e.g. yarp read) and press a button."   
4. At this moment launch gazebo (with the worldPROMPS.sdf if you have done movement from this world:   
gazebo -slibgazebo_yarp_clock.so worldPROMPS.sdf   
4.b.(new) use the command:  wholeBodyDynamicsTree --autoconnect --robot icubGazeboSim
5. Use the command:  iKinCartesianSolver --robot icubGazeboSim --part left_arm (from the path where the .ini linked to gazebo are)   
6. Use the command: simCartesianControl --robot icubGazeboSim (from the path where the .ini linked to gazebo are)   
7. Launch the programm replayTrajectories in ./CppProgram/build/bin   
8. Go back to the matlab windows and follow the instruction.   


## Actual development (V2)
In this version, we began to implement the code to be launch on the real iCub.   
We take into account information about external forces.   
We just need to be able to modify the compliance of the robot arm to finish this work.   
