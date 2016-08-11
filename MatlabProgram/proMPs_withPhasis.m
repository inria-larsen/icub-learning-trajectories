%In this programm we add to FourthAppliPromp the possibility to create
%trajectory with different velocity by computing the phasis.

close all;
clear all;

%variables
z = 100; %total time after modulation
nbFunctions(1) = 11; %number of basis functions
nbFunctions(2) = 21; %number of basis functions for forces
nbDof(1) = 3; %number of degree of freedom
nbDof(2) = 3; %number of forces
nbData = 30;%floor(2*z /3); %number of data before trying to find the correct movement

%Launch that only if you want to test it onto gazebo
%port open: port(/matlab/write)
%bottle b to write, c to read
initConnection;


%function test to create false data
%createVariables;

%we define also the variables nbKindOfTraj, var, totalTime
recoverTrajectories;


%compute the distribution for each kind of trajectories.
%we define var and TotalTime in this function
%here we need to define the bandwith of the gaussians h
h = 0.003; %bandwidth of the gaussians
computeDistributions;

%In this function we will play the trajectory into gazebo.
%TO do that launch 
%1. yarpserver
%2. gazebo with worldPROMPS
%3. go to a directory that contains the simCartesianLeftArm.ini /
%cartesianSolver.ini that have been modified to correspond to the gazebo
%simulation robot.
%4. Launch: simCartesianControl --robot icubGazeboSim in this terminal
%5. Create another terminal and launch:  iKinCartesianSolver --robot icubGazeboSim --part left_arm
%6. Then, launch the program. It will show you the learned distribution.
%7. Launch test_replay programm when the program is waiting for a
%connection.
%8 connect the port /matlab/write and /replay/read in the two sens.
%9 On matlab, give the number of the trajectory to replay and watch/
replay;


%creation of the begin of the movement.
%createMovementTest

%Recognition of the movement
recognitionTrajectory;

closeConnection;
