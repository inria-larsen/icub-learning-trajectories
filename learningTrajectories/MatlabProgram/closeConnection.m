%close the cpp program (that give order to gazebo)
disp('We close the module');
b.clear();
b.addDouble(-1);    
port.write(b);
port.close;
%portForces.close;
clear all;