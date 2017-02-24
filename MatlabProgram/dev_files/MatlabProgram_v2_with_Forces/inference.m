%in this function, we try to recongize a movement from the 30 first data
%and to complete it.
%TO do that, we consider the phasis of the movement as the mean of the
%phasis used during the learning.

%variable tuned to achieve the trajectory correctly
%accuracy = 0.000000001;
accuracy = 0.00000005;

%nbData = floor((mu_alpha(1) +mu_alpha(2) + mu_alpha(3))/0.03)
nbData = 60;
num= zeros(nbData,6);
%AskForData
b.clear();
b.addString('request_data');
b.addDouble(nbData);
port.write(b);
disp('Have send the message.');
c.clear();
port.read(c);
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

 y_trial = [num(1:nbData,1) ; num(1:nbData,2) ; num(1:nbData,3) ; filter(v,a,num(1:nbData,4)) ; filter(v,a,num(1:nbData,5)) ; filter(v,a,num(1:nbData,6))];
 y_trial_Tot = [num(:,1) ; num(:,2) ; num(:,3) ; filter(v,a,num(:,4)) ; filter(v,a,num(:,5)) ; filter(v,a,num(:,6))];

 
 affich = 1;
if(affich ==1) namePlot = ('x');
elseif(affich ==2) namePlot = ('y');
elseif(affich ==3) namePlot = ('z');
elseif(affich ==4) namePlot = ('fx');
elseif(affich ==5) namePlot = ('fy');
else namePlot = ('fz');
end
 % % Plot the trials we will use to verify the learning
nm = figure;
nm = visualisation(y_trial_Tot, 6, totalTimeTrial,affich, [1, 0, 0], nm);hold on;
nm = visualisation(y_trial, 6, nbData,affich, '.r', nm);
title('Trial trajectory');
xlabel('Iterations');
ylabel(namePlot);
%plot 3D trials
nam4 = figure;
nam4 = visualisation3D(y_trial_Tot, 6, totalTimeTrial(1), 0, [1, 0, 0], nam4);
nam4 = visualisation3D(y_trial, 3, nbData, 0, '.r', nam4);
title('Trial trajectory');
xlabel('x')
ylabel('y')
zlabel('z')
 

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
reco = {0 , -Inf, 1};

%% compute log_p(alpha) for a number of alphas
n_alpha_samples = 50;

%% defining deterministically a number of alphas between the minimum and the maximum observed
sampled_alphas = linspace(min_alpha, max_alpha, n_alpha_samples)';

%%
%log_p_alphas = log(mvnpdf(sampled_alphas,mu_alpha,Sigma_alpha));
%% mvnpdr retourne la densité de proba d'une distribution normal multivariable avec centre mu_alpha, et cov sigma_alpha

for k=1:50
    for i=1:nbKindOfTraj
        %matrix of cartesian basis functions that correspond to the first nbData 
        PSI_coor{k}{i} = computeBasisCoord(z,nbFunctions(1),sampled_alphas(k), floor(z/sampled_alphas(k)), h, nbData);
        %matrix of forces basis functions that correspond to the first nbData
        PSI_forces{k}{i} = computeBasisForces(z,nbFunctions(2),sampled_alphas(k), floor(z/sampled_alphas(k)), h, nbData);
        %matrix of basis functions for all data that correspond to the first
        %nbData
        PSI_mean{k}{i} = blkdiag(PSI_coor{k}{i},PSI_forces{k}{i});

        %we retrieve the learned trajectory of cartesian position
        u{k}{i} = PSI_coor{k}{i}*mu_w_coord{i};
        sigma{k}{i} = PSI_coor{k}{i}*sigma_w_coord{i}*PSI_coor{k}{i}' + accuracy*eye(size(PSI_coor{k}{i}*sigma_w_coord{i}*PSI_coor{k}{i}'));
        %we compute the probability it correspond to the actual trial
        for j=10:nbData
            prob{k}{i}(j)= my_log_mvnpdf([num(1:j,1) ; num(1:j,2) ; num(1:j,3)]',u{k}{i}(1:3*j)',sigma{k}{i}(1:3*j,1:3*j));
        end
      
        %we record the max of probability to know wich distribution we
        %recognize
        [ prob1, indx ] = max(prob{k}{i});
       
        if(prob1 > reco{2})
            reco{2} = prob1;
            reco{1} = i;
            reco{3} = k;
        end

    %     figure;
    %     visualisation(u{i}, 1, nbData, 'r');
    %     visualisation(u{i} + 1.96*sqrt(diag(sigma{i})), 1, nbData, ':r');
    %     visualisation(u{i} - 1.96*sqrt(diag(sigma{i})), 1, nbData, ':r');
    %     visualisation(y_trial{trial}, 1, nbData, 'b');
    end
end

quit = 0;
disp ''
disp(['The recognize trajectory is the number ', num2str(reco{1}), ' with p= ', num2str(reco{2}) ' and alpha= ', num2str(sampled_alphas(reco{3}))])
realAlpha = z /totalTimeTrial;
display(['The real alpha is ', num2str(realAlpha), ' with total time : ', num2str(totalTimeTrial) ])
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot 3D data
    nam2 = figure;
    for i=1:var(1)
        nam2 = visualisation3D(y{1}{i}, 6, totalTime(1,i), 0, [1, 0, 0], nam2);
    end
    for i=1:var(2)
        nam2 = visualisation3D(y{2}{i}, 6, totalTime(2,i), 0,[0, 0, 1], nam2);
    end
    for i=1:var(3)
        nam2 = visualisation3D(y{3}{i}, 6, totalTime(3,i),0,[0, 1, 0], nam2);
    end
    nam2 = visualisation3D(y_trial, 3, nbData, 0, 'v', nam2);
    nam2 = visualisation3D(y_trial_Tot, 3, totalTimeTrial, 0, ':', nam2);

    title('Trajectories used to learn')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    legend(nam2([2 (var(1)+2) (var(1) + var(2) + 2)]),'Right','Ahead','Top');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(reco{1} ==0) % if we don't recognizz the movement
    disp('We dont recognize the movement');
    msg = input('Send q to quit\n', 's');
    b.clear();
    b.addDouble(0.0)
    port.write(b);
    if(msg == 'q')
         disp('End of the programm.');  
    else
        clear mu_w_coord mu_w_f sigma_w_coord PSI_coor PSI_forces PSI_mean u sigma prob reco mu_n sigma_n y_trial realAlpha PSI_update ynew K mu_new sigma_new
        inference;
    end
else %if we recognize the movemet
    %we retrieve the computed distribution that correspond to the recognized
    %trajctory
    mu_new = mu_w{reco{1}};
    sigma_new = sigma_w{reco{1}}; 

    %we complete the data with the supposed forces correlated to the movement
    %according to the learned trajectory
    %y_trial =y_trial{trial}; %[y_trial{trial} ; PSI_forces{reco{1}}*mu_w_f{reco{1}}];

    %we suppose that the velocity of the movement is equal to the mean of
    %the one learned
    display(['alpha real: ', num2str(mu_alpha(reco{1})), ' totTimeTrial: ', num2str(z / mu_alpha(reco{1})) ])
    totalTimeExpexted = z / sampled_alphas(reco{3});
    display(['alpha_expected: ', num2str(sampled_alphas(reco{3})), ' totTime_exp: ', num2str(totalTimeExpexted) ])

    
    %we need to have the psi matrix and vector value according to time to
    %update the distribution (just a rewriting of data to simplify the next
    %computation.
    for t=1:nbData
        for i=1: nbDof(1) %+ nbDof(2)
            PSI_update{t}(i,:) = PSI_mean{reco{3}}{reco{1}}(t + nbData*(i-1),:);
            ynew{t}(i) = y_trial(t + nbData*(i-1)) ;
        end
    end
    
    % compute the new distribution (we try to pass by via point that correspond
    % to the fist nbData, with an accuracy tuned at the begining)
    for t=1:nbData     
        K = sigma_new*PSI_update{t}' * inv(accuracy*eye(size(PSI_update{t}*sigma_new*PSI_update{t}')) + PSI_update{t}*sigma_new*PSI_update{t}');
        mu_new = mu_new + K* (ynew{t}' - PSI_update{t}*mu_new);
        sigma_new = sigma_new - K*(PSI_update{t}*sigma_new);
    end

    clear PSI_update PSI_mean;
    b.clear();
    b.addString('send_data');
    port.write(b);
    
    
        % %%%%%%%%%%%%%%%%%%%%%%%%REPRESENTATION (plot)%%%%%%%%%%%%%%%
        %Plot the total trial and the data we have
    nameFig = figure;
    i = reco{1};

    for t=1:nbData
           nameFig(t) = scatter(t*realAlpha, y_trial((i-1)*nbData + t), '.b'); hold on;       
    end
    nameFig = visualisation2(y_trial_Tot, 6, totalTimeTrial, i, ':b', realAlpha, nameFig);

    %draw the infered movement
    nameFig = visualisation(PSI_z*mu_new, 6, z, i, '+g', nameFig);
    nameFig = visualisation(PSI_z*mu_w{i}, 6, z, i,'r', nameFig);
    nameFig = visualisation(PSI_z*(mu_w{i} + 1.96*sqrt(diag(sigma_w{i}))), 6, z,i, '-.r', nameFig);
    nameFig = visualisation(PSI_z*(mu_w{i}- 1.96*sqrt(diag(sigma_w{i}))), 6, z, i, '-.r', nameFig);
    legend(nameFig([1 (nbData+1) (nbData +2) (nbData +3) ]),'Data known', 'Data we should have', 'Data deducted', 'Learned distribution');
    trial = input('what movement have you tried? (1. left 2. ahead 3. top')
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
           name3D(t) = scatter3(y_trial(t), y_trial(nbData + t), y_trial(2*nbData + t), '.b'); hold on;
    end
    name3D = visualisation3D(y_trial_Tot, 6, totalTimeTrial,0, ':b', name3D);hold on;
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


    %draw the infered movement 3D forces
    name3D2 = figure;
    for t=1:nbData
           name3D2(t) = scatter3(y_trial(t + nbData*3), y_trial(nbData*4 + t), y_trial(5*nbData + t), '.b'); hold on;
    end
    %name3D2 = visualisation3D(y_trial_Tot{reco{1}}, 6, totalTimeExpexted,1, ':b', name3D2);hold on;
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

    replayRecognition;
end

