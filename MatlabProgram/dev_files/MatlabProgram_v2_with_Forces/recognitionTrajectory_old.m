%in this function, we try to recongize a movement from the 60 first data
%and to complete it.
%TO do that, we consider the phasis of the movement as the mean of the
%phasis used during the learning.
%we compute the new distribution without tacking account forces (that are
%C/C from the old distrib' to the new one.
%since we conserve in the sigma the corelation between the forces (where
%the data seems to be not normal) with the position, the standard deviation
%is not good
%variable tuned to achieve the trajectory correctly
accuracy = 0.00000001;

%accuracy that we want : choose randomly
trial = input('Give the test you want to do (1, 2, 3)\n');
disp(['we try the number ', num2str(trial)])

%computation of the loglikelihood for each trajectory using only cartesian
%coordinates

%we cut the mu_w to correspond only to the cartesian position informaiton
for i=1:nbKindOfTraj
    mu_w_coord{i} = mu_w{i}(1:nbDof(1)*nbFunctions(1));
    mu_w_f{i} = mu_w{i}(nbDof(1)*nbFunctions(1)+1:nbDof(1)*nbFunctions(1)+nbDof(2)*nbFunctions(2));
    sigma_w_coord{i} = sigma_w{i}(1:nbDof(1)*nbFunctions(1),1:nbDof(1)*nbFunctions(1));
end

% we compute for each learned distribution the loglikelihood that this
% movement correspond to the distribution
reco = {0 , -Inf };
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
    prob{i}= my_log_mvnpdf(y_trial{trial}(1:nbData*3)',u{i}',sigma{i});
    
    %we record the max of probability to know wich distribution we
    %recognize
    if(prob{i} > reco{2})
        reco{2} = prob{i};
        reco{1} = i;
    end
    
%     figure;
%     visualisation(u{i}, 1, nbData, 'r');
%     visualisation(u{i} + 1.96*sqrt(diag(sigma{i})), 1, nbData, ':r');
%     visualisation(u{i} - 1.96*sqrt(diag(sigma{i})), 1, nbData, ':r');
%     visualisation(y_trial{trial}, 1, nbData, 'b');
end

disp(['The recognize trajectory is the number ', num2str(reco{1})])
prob{1}
prob{2}
prob{3}
%we retrieve the computed distribution that correspond to the recognized
%trajctory
mu_n = mu_w{reco{1}}(1:nbDof(1)*nbFunctions(1));
sigma_n = sigma_w{reco{1}}(1:nbDof(1)*nbFunctions(1), (1:nbDof(1)*nbFunctions(1))); 

%we complete the data with the supposed forces correlated to the movement
%according to the learned trajectory
y_trial_nbData =y_trial{trial}; %[y_trial{trial} ; PSI_forces{reco{1}}*mu_w_f{reco{1}}];

%we aren't suppose to know "realData",  it is only used to draw the real
%trajectory of the sample if we continue it to the end
realAlpha = z /totalTimeTrial(trial);
display(['The real alpha is ', num2str(realAlpha), ' with total time : ', num2str(totalTimeTrial(reco{1})) ])
display(['The supposed alpha is ', num2str(mu_alpha(reco{1})), ' with total time : ', num2str(z / mu_alpha(reco{1})) ])


%we need to have the psi matrix and vector value according to time to
%update the distribution (just a rewriting of data to simplify the next
%computation.
for t=1:nbData
    for i=1: nbDof(1) %+ nbDof(2)
        PSI_update{t}(i,:) = PSI_mean{reco{1}}(t + nbData*(i-1),1:nbDof(1)*nbFunctions);
        ynew{t}(i) = y_trial_nbData(t + nbData*(i-1)) ;
    end
end

% compute the new distribution (we try to pass by via point that correspond
% to the fist nbData, with an accuracy tuned at the begining)
for t=1:nbData     
    K = sigma_n*PSI_update{t}' * inv(accuracy*eye(size(PSI_update{t}*sigma_n*PSI_update{t}')) + PSI_update{t}*sigma_n*PSI_update{t}');
    mu_n = mu_n + K* (ynew{t}' - PSI_update{t}*mu_n);
    sigma_n = sigma_n - K*(PSI_update{t}*sigma_n);
end
%we add old forces data
mu_new = [mu_n ; mu_w{reco{1}}(nbDof(1)*nbFunctions(1)+1:nbDof(1)*nbFunctions(1)+nbDof(2)*nbFunctions(2))];
sigma_new = [sigma_n, sigma_w{reco{1}}(1:nbDof(1)*nbFunctions(1) , nbDof(1)*nbFunctions(1)+1:nbDof(2)*nbFunctions(2)) ; sigma_w{reco{1}}(nbDof(1)*nbFunctions(1)+1:nbDof(1)*nbFunctions(1)+nbDof(2)*nbFunctions(2))];
    


%%%%%%%%%%%%%%%%%%%%%%%%REPRESENTATION (plot)%%%%%%%%%%%%%%%
%Plot the total trial and the data we have
nameFig = figure;
for t=1:nbData
    i=reco{1};
    %for i=1:nbDof(1) 
       nameFig(t) = scatter(t*realAlpha, y_trial_nbData((i-1+3)*nbData + t), '.b'); hold on;       
    %end
end
nameFig = visualisation2(y_trial_Tot{reco{1}}, 6, totalTimeTrial(reco{1}), i+3, ':b', realAlpha, nameFig);
%draw the infered movement

nameFig = visualisation(PSI_z*mu_new, 6, z,i+3, '+g', nameFig);
nameFig = visualisation(PSI_z*mu_w{i}, 6, z,i+3, 'r', nameFig);
nameFig = visualisation(PSI_z*(mu_w{i} + 1.96*sqrt(diag(sigma_w{i}))), 6, z, i+3,'-.r', nameFig);
nameFig = visualisation(PSI_z*(mu_w{i}- 1.96*sqrt(diag(sigma_w{i}))), 6, z, i+3, '-.r', nameFig);
legend(nameFig([1 (nbData+1) (nbData +2) (nbData +3) ]),'Data known', 'Data we should have', 'Data deducted', 'Learned distribution');
if(trial ==1)
    name = 'right';
elseif(trial == 2)
    name = 'ahead';
else
    name= 'top';
end
title(['Position recognition of the ', name ,' trajectory ']);
xlabel('Iterations');
ylabel('x cartesian position (m)');

%draw the infered movement 3D
name3D = figure;
for t=1:nbData
       name3D(t) = scatter3(y_trial_nbData(t), y_trial_nbData(nbData + t), y_trial_nbData(2*nbData + t), '.b'); hold on;
end
name3D = visualisation3D(y_trial_Tot{reco{1}}, 6, totalTimeTrial(reco{1}),0, ':b', name3D);hold on;
%nam4 = visualisation3D(y_trial_Tot{1}, 6, totalTimeTrial(1), 0, [1, 0, 0], nam4);


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
clear mu_w_coord mu_w_f sigma_w_coord PSI_coor PSI_forces PSI_mean u sigma prob mu_n sigma_n y_trial_nbData realAlpha PSI_update ynew K mu_new sigma_new reco

%replayRecognition;