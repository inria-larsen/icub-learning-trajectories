%in this function, we try to recongize a movement from the 30 first data
%and to complete it.
%TO do that, we consider the phasis of the movement as the mean of the
%phasis used during the learning.

%variable tuned to achieve the trajectory correctly

%accuracy that we want : choose randomly
% text= ['Give the test you want to do (from 1 to ' num2str(nbKindOfTraj) ')'];
% trial = input(text);
% disp(['we try the number ', num2str(trial)])
 accuracy=0.00000000000000001;
%while(accuracy <= 0.1)
 %pb ici
%y_trial_Tot{1} = ymean{1};
%totalTimeTrial(1) =100;% totalTime(1,5);
trial=1;

%computation of the loglikelihood for each trajectory using only cartesian
%coordinates
typeReco = 1; %what kind of dof we use for the reco{cpt}
vr = 0;
prevDof = 0;
for j=1:typeReco-1
    vr = vr + nbDof(j)*nbFunctions(j);
    prevDof = prevDof + nbDof(j); 

end
%we cut the mu_w to correspond only to the cartesian position informaiton
for i=1:nbKindOfTraj
    mu_w_reco{i} = mu_w{i}(vr + 1 : vr + nbDof(typeReco)*nbFunctions(typeReco));
    sigma_w_reco{i} = sigma_w{i}(vr + 1 : nbDof(typeReco)*nbFunctions(typeReco),vr + 1:nbDof(typeReco)*nbFunctions(typeReco));
end


% we compute for each learned distribution the loglikelihood that this
% movement correspond to the distribution

%reco = cell(11,3);
meanTime= floor((z / mu_alpha));
cpt = 1;



%psiTrial = computeBasisFunction(z,nbFunctions, nbDof, mu_alpha(1), floor(z/mu_alpha(1)), center_gaussian, h,totalTimeTrial(trial));
%for step=1:floor(meanTimeDiv2/10):meanTimeDiv2

for step=1: floor(totalTimeTrial(1) /10): totalTimeTrial(1)
    reco{cpt} = [0.0 , 0.0, 0.0 ];
    for i=1:nbKindOfTraj
        cpt
        %matrix of basis functions that correspond to the first nbData for
        %the DOF taked into account for the recognition
        PSI_reco{step}{i} = computeBasisFunction(z,nbFunctions(typeReco), nbDof(typeReco), mu_alpha(i), floor(z/mu_alpha(i)), center_gaussian(typeReco), h(typeReco),step);
        %matrix of basis functions for all data that correspond to the first
        %nbData
        PSI_tot{step}{i} = computeBasisFunction(z,nbFunctions, nbDof, mu_alpha(i), floor(z/mu_alpha(i)), center_gaussian, h,step);

        %we retrieve the learned trajectory of cartesian position
        u{step}{i} = PSI_reco{step}{i}*mu_w_reco{i};
        sigma{step}{i} = PSI_reco{step}{i}*sigma_w_reco{i}*PSI_reco{step}{i}' + accuracy*eye(size(PSI_reco{step}{i}*sigma_w_reco{i}*PSI_reco{step}{i}'));
        
        %we compute the probability it correspond to the actual trial
        y_trial_part{step} = [];
        for k=1:nbDof(typeReco)
           y_trial_part{step} = [y_trial_part{step} ; y_trial_Tot{trial}(totalTimeTrial(trial)*(prevDof+k -1) + 1 : totalTimeTrial(trial)*(prevDof+k -1) + step)];
        end
        prob(i,cpt)= logLikelihood(y_trial_part{step}',u{step}{i}',sigma{step}{i})

        %we record the max of probability to know wich distribution we
        %recognize
    %    if(prob(i,step) > reco{cpt}(2))
            reco{cpt}(2) = prob(i,cpt);
            reco{cpt}(1) = i;
            reco{cpt}(3) = step;
     %   end

    end

    disp(['The recognize trajectory is the number ', num2str(reco{cpt}(1))])


    %we retrieve the computed distribution that correspond to the recognized
    %trajctory
    mu_new = mu_w{reco{cpt}(1)};
    sigma_new = sigma_w{reco{cpt}(1)}; 
    %logLikelihood( y_trial_Tot{trial}', mu_new', sigma_new) ;
    %we complete the data with the supposed forces correlated to the movement
    %according to the learned trajectory


    y_prev = [];
    for i=1:prevDof
        y_prev  = [y_prev ; y_trial_Tot{trial}(totalTimeTrial(trial)*(i -1) + 1 : totalTimeTrial(trial)*(i -1) + reco{cpt}(3))];
    end

    y_after = [];
    for i= (1 + prevDof + nbDof(typeReco)): nbDofTot
       y_after = [y_after;  y_trial_Tot{trial}(totalTimeTrial(trial)*(i -1) + 1 : totalTimeTrial(trial)*(i -1) + reco{cpt}(3))];   
    end
    %all data (all DOF) from t=1 to the time of inference
    y_trial_nbData =[y_prev ; y_trial_part{reco{cpt}(3)}; y_after ];

    %we aren't suppose to know "realData",  it is only used to draw the real
    %trajectory of the sample if we continue it to the end
    realAlpha = z /totalTimeTrial(reco{cpt}(1));
 %   display(['The real alpha2 is ', num2str(realAlpha), ' with total time : ', num2str(totalTimeTrial(reco{cpt}(1))) ])
 %   display(['The supposed alpha2 is ', num2str(mu_alpha(reco{cpt}(1))), ' with total time : ', num2str(z / mu_alpha(reco{cpt}(1))) ])


    %we need to have the psi matrix and vector value according to time to
    %update the distribution (just a rewriting of data to simplify the next
    %computation.
    vr=0;
    for j=1:typeReco-1
        vr = vr + nbDof(j);
    end
    for t=1:reco{cpt}(3)
        for i=vr + 1: vr + nbDof(typeReco)
            PSI_update{t}(i,:) = PSI_tot{reco{cpt}(3)}{reco{cpt}(1)}(t + reco{cpt}(3)*(i-1),:);
            ynew{t}(i) = y_trial_nbData(t + reco{cpt}(3)*(i-1)) ;
        end
    end

%fig = figure(1)

%filename = ['test', num2str(cpt), '.gif'];
    % compute the new distribution (we try to pass by via point that correspond
    % to the fist nbData, with an accuracy tuned at the begining)
    for t=1:reco{cpt}(3)     
        K = sigma_new*PSI_update{t}' * inv(accuracy*eye(size(PSI_update{t}*sigma_new*PSI_update{t}')) + PSI_update{t}*sigma_new*PSI_update{t}');
        mu_new = mu_new + K* (ynew{t}' - PSI_update{t}*mu_new);
        sigma_new = sigma_new - K*(PSI_update{t}*sigma_new);
%        if(step == totalTimeTrial(1)) 
%        
%         shadedErrorBar([1:1:100],PSI_z*mu_new, PSI_z*1.96*sqrt(diag(sigma_new)), 'g');
%         end
%         fig(size(fig,2)+1) = plot([1:1:100],PSI_z*mu_new,'g');
%         fig =  visualisation2(y_trial_Tot{i}, nbDofTot, totalTimeTrial(i),1, '--magenta', realAlpha, fig);hold on;
% 
%          %visualisation2(y_trial_Tot{i}, nbDofTot, totalTimeTrial(i),type, '-.r', mu_alpha(i), nf3D);hold on;
%          %legend(fig([(size(fig,2)-1) size(fig,2)]),'inference', 'real trajectory')
%          drawnow
%          
%         frame = getframe(1);
% 
%         im = frame2im(frame);
%          [imind,cm] = rgb2ind(im,256);
%          if( t == 1)
%               imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%          else
%              imwrite(imind,cm,filename,'gif','WriteMode','append');
%          end
       
         %input('wait');
    end

    newmu{cpt} = mu_new;
    newSigma{cpt} = sigma_new;
    cpt = cpt+1;
    %reco{cpt} = [reco{cpt-1}(1),reco{cpt-1}(2),reco{cpt-1}(3)]; %{trajectory recognized, probability, timestep} of the best likelihood
%     
%     clear u ;
%     u = psiTrial*mu_new;
%     S = psiTrial*sigma_new*psiTrial' + accuracy*eye(size(psiTrial*sigma_new*psiTrial'));
%     logkike(step) = logLikelihood( y_trial_Tot{trial}', u', S); 
    clear mu_new sigma_new u psiTrial;
end
    Drawing;
%accuracy = accuracy*10;
%end
%clear cpt mu_w_coord mu_w_f sigma_w_coord PSI_coor PSI_forces PSI_mean u sigma prob reco mu_n sigma_n y_trial_nbData realAlpha PSI_update ynew K mu_new sigma_new