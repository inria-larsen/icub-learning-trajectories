LoadYarp;
import yarp.Port;
import yarp.Bottle;

done=0;

port=Port;
%first close the port just in case
port.close;

disp('Going to open port /matlab/write and read');
port.open('/matlab/write');
rep = input('Please connect to a bottle sink (e.g. yarp read) and press a button.\n');
b = Bottle;
c = Bottle;