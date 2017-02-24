%This toolbox allow (1) to learn the distribution of any kind of trajectories (i.e. data input that evolved in time) (2) to infer the end of an initiate trajectory, thanks to the learned distribution.
%by Oriane Dermy 07/09/2016
% For any problem / remark / improvement, contact me:
% oriane.dermy@gmail.com with subject [proMPs_toolbox]
close all;
clearvars;

%%%VARIABLES, please refer you to the readme
nbKindOfTraj = 3;
 nameD{1} = 'Data/dataAhead.txt';
 nameD{2} = 'Data/dataTop.txt';
 nameD{3} = 'Data/dataRight.txt';
z = 100; %total time after modulation

nbDof(1) = 3; %number of degree of freedom
nbDof(2) = 3; %number of forces

nbFunctions(1) = 5;%5 %51; %number of basis functions
nbFunctions(2) = 5; %21; %number of basis functions for forces
nbTotFunctions = 0;

for i=1:size(nbFunctions,2)
    nbTotFunctions = nbTotFunctions + nbFunctions(i)*nbDof(i);
end
center_gaussian(1) = 1.0 / (nbFunctions(1));
center_gaussian(2) = 1.0 / (nbFunctions(2));
h(1) = center_gaussian(1)/5;%0.02%center_gaussian(1)*(1/z); %0.006; %bandwidth of the gaussians
h(2) = center_gaussian(1)/5;%center_gaussian(1)*6*(1/z)/100;%0.003;
accuracy =0.00000005; %precision we want during inference


%Launch that only if you want to test it onto gazebo
%port open: port(/matlab/write)
%bottle b to write, c to read
initConnection;
    
    
nbData = 40 %floor(2*z /3); %number of data max with what you try to find the correct movement
%createToyExample

%we define also the variables nbKindOfTraj, var, totalTime
recoverData;

%plot recoverData
%draw_RecoverData

%compute the distribution for each kind of trajectories.
%we define var and TotalTime in this function
%here we need to define the bandwith of the gaussians h
%computeDistributions_withCrossOver;
computeDistributions

%plot distribution
%drawDistribution


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

%Recognition of the movement
inferenceFromZero;
%draw the infered movement
drawInferedMovement

replayRecognition;

%drawnow();
closeConnection;



