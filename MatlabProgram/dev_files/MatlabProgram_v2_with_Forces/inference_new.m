accuracy = 0.00000005;

%nbData = floor((mu_alpha(1) +mu_alpha(2) + mu_alpha(3))/0.03)
nbData = 60;

%AskForData
b.clear();
b.addString('request_data');
b.addDouble(nbData);
port.write(b);
disp('Have send the message.');
c.clear();
port.read(c);
disp('we have received data');
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

 y_trial_nbData = [num(1:nbData,1) ; num(1:nbData,2) ; num(1:nbData,3) ; filter(v,a,num(1:nbData,4)) ; filter(v,a,num(1:nbData,5)) ; filter(v,a,num(1:nbData,6))];
 y_trial_Tot = [num(:,1) ; num(:,2) ; num(:,3) ; filter(v,a,num(:,4)) ; filter(v,a,num(:,5)) ; filter(v,a,num(:,6))];

 %%%%%%%%%%%%%DEBUG
   % Plot the learned distribution for the x cartesian position of the
   % trajectory i
   %nf = visualisation(PSI_z*mu_w{1}, 6, z, 1, 'r', nf);
   %nf = visualisation(PSI_z*(mu_w{1} + 1.96*sqrt(diag(sigma_w{1}))), 6, z, 1,'-.r', nf);
   %nf = visualisation(PSI_z*(mu_w{1}- 1.96*sqrt(diag(sigma_w{1}))), 6, z, 1, '-.r', nf);
   for varr=1:3
      nf = figure;         
      for j=1:var(varr)
       nf = visualisation2(y{varr}{j}, 6, totalTime(varr,j),varr, 'b', alpha{varr}(j) , nf);hold on;
      end 
         nf = visualisation2(y_trial_Tot, 6, totalTimeTrial,varr, 'k', 100/totalTimeTrial, nf);hold on;
    
     title('learned distribution for the x position of the trajectory right');%strcat('learned distribution for the', atostr(i), ' trajectory');
    
    xlabel('Iteration')
    ylabel('x (m)')
   legend(nf([2 4 5 size(nf,2)]),'Distribution Learned','Standard deviation','Data', 'trial');   
   end

 %%%%%%%%%%%%%%%%%END DEBUG 
 for i=1:nbKindOfTraj
    mu_w_coord{i} = mu_w{i}(1:nbDof(1)*nbFunctions(1));
    mu_w_f{i} = mu_w{i}(nbDof(1)*nbFunctions(1)+1:nbDof(1)*nbFunctions(1)+nbDof(2)*nbFunctions(2));
    sigma_w_coord{i} = sigma_w{i}(1:nbDof(1)*nbFunctions(1),1:nbDof(1)*nbFunctions(1));
end

% we compute for each learned distribution the loglikelihood that this
% movement correspond to the distribution
reco = {0, -Inf };

figure;
for i=1:nbKindOfTraj

    %matrix of cartesian basis functions that correspond to the first nbData 
    PSI_coor{i} = computeBasisCoord(z,nbFunctions(1),mu_alpha(i), floor(z/mu_alpha(i)), h, nbData);
    %matrix of forces basis functions that correspond to the first nbData
    PSI_forces{i} = computeBasisForces(z,nbFunctions(2),mu_alpha(i), floor(z/mu_alpha(i)), h, nbData);
    %matrix of basis functions for all data that correspond to the first
    %nbData
    PSI_mean{i} = blkdiag(PSI_coor{i},PSI_forces{i});
    
    %we retrieve the learned trajectory of cartesian position
    u{i} = PSI_coor{i}*mu_w_coord{i};
    sigma{i} = PSI_coor{i}*sigma_w_coord{i}*PSI_coor{i}' + accuracy*eye(size(PSI_coor{i}*sigma_w_coord{i}*PSI_coor{i}'));
    %we compute the probability it correspond to the actual trial
    %TODO changer ça
    prob{i}= -mean(abs(y_trial_nbData(1:nbData*3) -u{i})) %my_log_mvnpdf(y_trial_nbData(1:nbData*3)',u{i}',sigma{i});
    plot(abs(y_trial_nbData(1:nbData*3) -u{i}), 'k');hold on;
    plot(y_trial_nbData(1:nbData*3),'b');
    if(i==1) plot(u{i},'r');
    elseif(i==2) plot(u{i},'g');
    else plot(u{i},'m');
    end
   
    %we record the max of probability to know wich distribution we
    %recognize
    if(prob{i} > reco{2})
        reco{2} = prob{i};
        reco{1} = i;
    end
 
end

disp(['The recognize trajectory is the number ', num2str(reco{1})])
prob{1}
prob{2}
prob{3}
%we retrieve the computed distribution that correspond to the recognized
%trajctory
mu_new = mu_w{reco{1}};
sigma_new = sigma_w{reco{1}}; 

%we aren't suppose to know "realData",  it is only used to draw the real
%trajectory of the sample if we continue it to the end
realAlpha = z /totalTimeTrial;
display(['The real alpha is ', num2str(realAlpha), ' with total time : ', num2str(totalTimeTrial) ])
display(['The supposed alpha is ', num2str(mu_alpha(reco{1})), ' with total time : ', num2str(z / mu_alpha(reco{1})) ])


%we need to have the psi matrix and vector value according to time to
%update the distribution (just a rewriting of data to simplify the next
%computation.
for t=1:nbData
    for i=1: nbDof(1) %+ nbDof(2)
        PSI_update{t}(i,:) = PSI_mean{reco{1}}(t + nbData*(i-1),:);
        ynew{t}(i) = y_trial_nbData(t + nbData*(i-1)) ;
    end
end

% compute the new distribution (we try to pass by via point that correspond
% to the fist nbData, with an accuracy tuned at the begining)
for t=1:nbData     
    K = sigma_new*PSI_update{t}' * inv(accuracy*eye(size(PSI_update{t}*sigma_new*PSI_update{t}')) + PSI_update{t}*sigma_new*PSI_update{t}');
    mu_new = mu_new + K* (ynew{t}' - PSI_update{t}*mu_new);
    sigma_new = sigma_new - K*(PSI_update{t}*sigma_new);
end

%%%%%%%%%%%%%%%%%%%%%%%%REPRESENTATION (plot)%%%%%%%%%%%%%%%
%Plot the total trial and the data we have
nameFig = figure;
i = reco{1};

for t=1:nbData
       nameFig(t) = scatter(t*realAlpha, y_trial_nbData((i-1)*nbData + t), '.b'); hold on;       
end
nameFig = visualisation2(y_trial_Tot, 6, totalTimeTrial, i, ':b', realAlpha, nameFig);
%draw the infered movement
nameFig = visualisation(PSI_z*mu_new, 6, z, i, '+g', nameFig);
nameFig = visualisation(PSI_z*mu_w{i}, 6, z, i,'r', nameFig);
nameFig = visualisation(PSI_z*(mu_w{i} + 1.96*sqrt(diag(sigma_w{i}))), 6, z,i, '-.r', nameFig);
nameFig = visualisation(PSI_z*(mu_w{i}- 1.96*sqrt(diag(sigma_w{i}))), 6, z, i, '-.r', nameFig);
legend(nameFig([1 (nbData+1) (nbData +2) (nbData +3) ]),'Data known', 'Data we should have', 'Data deducted', 'Learned distribution');
trial = input('what do you tried ? 1. right 2. ahead 3. top');
if(trial ==1)
    name = 'right';
    name2 = 'y cartesian position' ;

elseif(trial == 2)
    name = 'ahead';
    name2 = 'x cartesian position';
else
    name= 'top';
    name2 = 'z cartesian position';
end
title(['Position recognition of the ', name ,' trajectory ']);
xlabel('Iterations');
ylabel(name2);

%draw the infered movement 3D
name3D = figure;
for t=1:nbData
       name3D(t) = scatter3(y_trial_nbData(t), y_trial_nbData(nbData + t), y_trial_nbData(2*nbData + t), '.b'); hold on;
end
name3D = visualisation3D(y_trial_Tot, 6, totalTimeTrial,0, ':b', name3D);hold on;


i = reco{1};
name3D = visualisation3D(PSI_z*mu_new, 6, z, 0, '+g', name3D);
name3D = visualisation3D(PSI_z*mu_w{i}, 6, z, 0, '.r', name3D);
name3D = visualisation3D(PSI_z*(mu_w{i} + 1.96*sqrt(diag(sigma_w{i}))), 6, z, 0, '-.r', name3D);
name3D = visualisation3D(PSI_z*(mu_w{i}- 1.96*sqrt(diag(sigma_w{i}))), 6, z, 0, '-.r', name3D);
legend(name3D([1 (nbData+1) (nbData +2) (nbData +3) ]),'Data known', 'Data we should have', 'Data deducted', 'Learned distribution');

if(trial ==1)
    name = 'right';
elseif(trial == 2)
    name = 'ahead';
else
    name= 'top';
end
title(['Position recognition of the ', name ,' trajectory'])
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');


%draw the infered movement 3D forces
name3D2 = figure;
for t=1:nbData
       name3D2(t) = scatter3(y_trial_nbData(t + nbData*3), y_trial_nbData(nbData*4 + t), y_trial_nbData(5*nbData + t), '.b'); hold on;
end
name3D2 = visualisation3D(y_trial_Tot, 6, totalTimeTrial,1, ':b', name3D2);hold on;
%nam4 = visualisation3D(y_trial_Tot{1}, 6, totalTimeTrial(1), 0, [1, 0, 0], nam4);


i = reco{1};
name3D2 = visualisation3D(PSI_z*mu_new, 6, z, 0, '+g', name3D2);
name3D2 = visualisation3D(PSI_z*mu_w{i}, 6, z, 0, '.r', name3D2);
name3D2 = visualisation3D(PSI_z*(mu_w{i} + 1.96*sqrt(diag(sigma_w{i}))), 6, z, 1, '-.r', name3D2);
name3D2 = visualisation3D(PSI_z*(mu_w{i}- 1.96*sqrt(diag(sigma_w{i}))), 6, z, 1, '-.r', name3D2);
legend(name3D2([1 (nbData+1) (nbData +2) (nbData +3) ]),'Data known', 'Data we should have', 'Data deducted', 'Learned distribution');

if(trial ==1)
    name = 'right';
elseif(trial == 2)
    name = 'ahead';
else
    name= 'top';
end
title(['Forces recognition of the ', name ,' trajectory'])
xlabel('fx (m)');
ylabel('fy (m)');
zlabel('fz (m)');
 

clear PSI_update PSI_mean;
    b.clear();
    b.addString('send_data');
    port.write(b);
    
replayRecognition
%clear mu_w_coord mu_w_f sigma_w_coord PSI_coor PSI_forces PSI_mean u sigma prob reco mu_n sigma_n y_trial_nbData realAlpha PSI_update ynew K mu_new sigma_new
