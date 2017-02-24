%AskForData
initConnection;
b.clear();
b.addString('request_data');
b.addDouble(30);
port.write(b);
disp('Have send the message.');
c.clear();
port.read(c);
disp('Have receive data.');
disp(c);

num2 = str2num(c);
totalTimeTrial =(size(num2,2)/6);
a = 1;
v = [1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10];
for(t=1:totalTimeTrial)
    num(t,1) = num2(6*(t-1) + 1);
    num(t,2) = num2(6*(t-1) + 2);
    num(t,3) = num2(6*(t-1) + 3);
    num(t,4) = num2(6*(t-1) + 4);
    num(t,5) = num2(6*(t-1) + 5);
    num(t,6) = num2(6*(t-1) + 6);
    %disp(['Receiving: x = ', num2str(num(t,1)), ', y = ',num2str(num(t,2)), ', z = ', num2str(num(t,3)), 'fx = ', num2str(num(t,4)), ', fy = ',num2str(num(t,5)), ', fz = ', num2str(num(t,6)) ]);
end


clear Bt;
load('dataTopTotal.mat', 'Bt')
Bt{size(Bt,2) +1} = num;
%Bt{1} = num;
save('dataTopTotal.mat', 'Bt')
closeConnection;
clear all;