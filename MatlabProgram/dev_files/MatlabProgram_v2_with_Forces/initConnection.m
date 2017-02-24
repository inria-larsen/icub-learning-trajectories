LoadYarp;
import yarp.Port;
import yarp.Bottle;

port=Port;
%first close the port just in case
port.close;
% portForces.close;
disp('Going to open port /matlab/write and read');
port.open('/matlab/write');

% portForces.open('/matlab/read');
%port.setTimeout(1);
rep = input('Please create a port (e.g. yarp read) and press a button.\n');
%disp 'Please create a port (e.g. yarp read) and press a button.\n'
disp('connecting matlab to verify')
while (yarp.Network.connect('/matlab/write','/verify/read', 'tcp')~=1)
    disp('...');
end
disp '[success] port ';

% if yarp.Network.connect('/wholeBodyDynamics/left_arm/ext_ft_sens:o','/matlab/read')  % autoconnect just for testing
%     disp '[success] port /wholeBodyDynamics/left_arm/ext_ft_sens:o connected to /matlab/read';
% else
%     disp '[warning] port NOT connected to /wholeBodyDynamics/left_arm/ext_ft_sens:o, does it exist?';
% end

%rep = input('Please connect to a bottle sink (e.g. yarp read) and press a button.\n');
b = Bottle;
c = Bottle;
% readForces = Bottle;
