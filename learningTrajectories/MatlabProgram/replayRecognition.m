%In this function we will play the trajectory recognized into gazebo.

data = PSI_z*mu_new;
data_max = PSI_z*(mu_new + 1.96*sqrt(diag(sigma_new)));


for t = round(mu_alpha(reco{1})*nbData): z
    b.clear();
    for i = 1 : nbDof(1) %+ nbDof(2)
        val(t,i) = data(z*(i-1)+t);
        b.addDouble(val(t,i));
    end
    dist_max = abs(data(z*(4-1)+t) - data_max(z*(4-1)+t)) + abs(data(z*(5-1)+t) - data_max(z*(5-1)+t)) + abs(data(z*(6 -1)+t) - data_max(z*(6 -1)+t));
    % message = num2string(val(t,:))
    port.write(b);
    %disp('Have send the message.');
    port.read(c);
    disp(c);
    num = str2num(c);
    disp(['Receiving: fx = ', num2str(num(1,1)), ', fy = ',num2str(num(1,2)), ', fy = ', num2str(num(1,3))]);
    disp(['Expected:  fx = ', num2str(data(z*(4-1)+t)), ', fy = ',num2str(data(z*(5-1)+t)), ', fz = ', num2str(data(z*(6 -1)+t))]);
    dist = abs(data(z*(4-1)+t) - num(1,1)) + abs(data(z*(5-1)+t) - num(1,2)) + abs(data(z*(6 -1)+t) - num(1,3));

    if(dist > dist_max)
        disp('The robot receive a too big force, so it has to become compliant!')
        dist 
        dist_max
        compliant = 0.0;
    else
        compliant = 1.0 - (dist / dist_max);
    end
    disp(['compliance = ', num2str(compliant)]);
end
msg = input('Send q to quit\n', 's');
if(msg == 'q')
    disp('End of the programm.');
else
    b.clear();
    b.addDouble(0.0);
    port.write(b);
    recognitionTrajectory;
end
