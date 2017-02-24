%in this function, we try to recongize a movement from the 30 first data
%and to complete it.
%TO do that, we consider the phasis of the movement as the mean of the
%phasis used during the learning.
%close all

   trial=1
%for trial=1:10
    %variable tuned to achieve the trajectory correctly
  %  clearvars -except max_alpha stock_tps trajErr trial stock_err newmu newSigma trajInf logTest  PSI_inf alpha_inf
  %  load('Data/dataDistributionManyDOF2.mat');
    
    %%
%     clear y_trial_Tot totalTimeTrial
%     y_trial_Tot =  y;
%     totalTimeTrial = totalTime;
    
    
    %accuracy that we want : choose randomly
    accuracy=0.0000000001;


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
    nbStep= 10; 
    for valstep=1: floor(totalTimeTrial(1,trial) /nbStep): min(totalTimeTrial(1, trial), z/ max_alpha)
        totstep = size([1: floor(totalTimeTrial(1,trial) /nbStep): min(totalTimeTrial(1, trial), z/ max_alpha)],2);
        tstart = tic;
        reco{cpt} = [0.0 , 0.0, 0.0 ];
        for i=1:nbKindOfTraj
            %we compute the probability it correspond to the actual trial of
            %"step" steps.
            y_trial_part{valstep} = [];
            for k=1:nbDof(typeReco)
               y_trial_part{valstep} = [y_trial_part{valstep} ; y_trial_Tot{1}{trial}(totalTimeTrial(1, trial)*(prevDof+k -1) + 1 : totalTimeTrial(1,trial)*(prevDof+k -1) + valstep)];
            end

            PSI_recoMin{valstep}{i} = computeBasisFunction(z,nbFunctions(typeReco), nbDof(typeReco), min_alpha(i), floor(z/min_alpha(i)), center_gaussian(typeReco), h(typeReco),valstep);
            PSI_recoMean{valstep}{i} = computeBasisFunction(z,nbFunctions(typeReco), nbDof(typeReco), mu_alpha(i), floor(z/mu_alpha(i)), center_gaussian(typeReco), h(typeReco),valstep);
            PSI_recoMax{valstep}{i} = computeBasisFunction(z,nbFunctions(typeReco), nbDof(typeReco), max_alpha(i), floor(z/max_alpha(i)), center_gaussian(typeReco), h(typeReco),valstep);
            distMin = mean(abs(PSI_recoMin{valstep}{i}*mu_w_reco{i} -y_trial_part{valstep}));
            distMean = mean(abs(PSI_recoMean{valstep}{i}*mu_w_reco{i} -y_trial_part{valstep}));
            distMax = mean(abs(PSI_recoMax{valstep}{i}*mu_w_reco{i} -y_trial_part{valstep}));
            sumDist =  distMin + distMax + distMean;
            wtestMin = (PSI_recoMin{valstep}{i}'*PSI_recoMin{valstep}{i}+ 1e-6*eye(size( PSI_recoMin{valstep}{i}'*PSI_recoMin{valstep}{i},1))) \ PSI_recoMin{valstep}{i}'* y_trial_part{valstep};    
            wtestMax = (PSI_recoMax{valstep}{i}'*PSI_recoMax{valstep}{i}+ 1e-6*eye(size( PSI_recoMax{valstep}{i}'*PSI_recoMax{valstep}{i},1))) \ PSI_recoMax{valstep}{i}'* y_trial_part{valstep};    
            wtestMean = (PSI_recoMean{valstep}{i}'*PSI_recoMean{valstep}{i}+ 1e-6*eye(size( PSI_recoMean{valstep}{i}'*PSI_recoMean{valstep}{i},1))) \ PSI_recoMean{valstep}{i}'* y_trial_part{valstep};    
           % la = logLikelihood(wtestMin,mu_w_reco{i},sigma_w_reco{i} );
           % lb = logLikelihood(wtestMean,mu_w_reco{i},sigma_w_reco{i} );
           %M lc = logLikelihood(wtestMax,mu_w_reco{i},sigma_w_reco{i} );
            %disp(['a: ', num2str(la), ' b: ', num2str(lb), ' c: ', num2str(lc) ]);

            alpha_inf{trial}(cpt) = ((distMax+ distMean) / (2*sumDist))*min_alpha(i) + ((distMin+ distMean) / (2*sumDist))*max_alpha(i) + ((distMax+ distMin) / (2*sumDist))*mu_alpha(i);
           % disp(['We supposed from alpha(min:mean:max) = [' num2str(min_alpha(i)), ':', num2str(mu_alpha(i)), ':', num2str(max_alpha(i)), '], that alpha = ', num2str(alpha_inf), '. Normaly should be: ', num2str(z/totalTimeTrial(trial)) ]);
            PSI_reco{valstep}{i} = computeBasisFunction(z,nbFunctions(typeReco), nbDof(typeReco), alpha_inf{trial}(cpt), floor(z/alpha_inf{trial}(cpt)), center_gaussian(typeReco), h(typeReco),valstep);

            PSI_tot{valstep}{i} = computeBasisFunction(z,nbFunctions, nbDof, alpha_inf{trial}(cpt), floor(z/alpha_inf{trial}(cpt)), center_gaussian, h,valstep);
            wtest = (PSI_reco{valstep}{i}'*PSI_reco{valstep}{i}+ 1e-6*eye(size( PSI_reco{valstep}{i}'*PSI_reco{valstep}{i},1))) \ PSI_reco{valstep}{i}'* y_trial_part{valstep};    
            prob(i,cpt) = logLikelihood(wtest,mu_w_reco{i},sigma_w_reco{i} );

            %we record the max of probability to know wich distribution we
            %recognize
        %    if(prob(i,step) > reco{cpt}(2))
                reco{cpt}(2) = prob(i,cpt);
                reco{cpt}(1) = i;
                reco{cpt}(3) = valstep;
         %   end

        end

        %disp(['The recognize trajectory is the number ', num2str(reco{cpt}(1))])


        %we retrieve the computed distribution that correspond to the recognized
        %trajctory
        mu_new = mu_w{reco{cpt}(1)};
        sigma_new = sigma_w{reco{cpt}(1)}; 
        %logLikelihood( y_trial_Tot{trial}', mu_new', sigma_new) ;
        %we complete the data with the supposed forces correlated to the movement
        %according to the learned trajectory


        y_prev = [];
        for i=1:prevDof
            y_prev  = [y_prev ; y_trial_Tot{1}{trial}(totalTimeTrial(1, trial)*(i -1) + 1 : totalTimeTrial(1,trial)*(i -1) + reco{cpt}(3))];
        end

        y_after = [];
        for i= (1 + prevDof + nbDof(typeReco)): nbDofTot
           y_after = [y_after;  y_trial_Tot{1}{trial}(totalTimeTrial(1, trial)*(i -1) + 1 : totalTimeTrial(1, trial)*(i -1) + reco{cpt}(3))];   
        end
        %all data (all DOF) from t=1 to the time of inference
        y_trial_nbData =[y_prev ; y_trial_part{reco{cpt}(3)}; y_after ];

        %we aren't suppose to know "realData",  it is only used to draw the real
        %trajectory of the sample if we continue it to the end
        realAlpha = z /totalTimeTrial(1, reco{cpt}(1));
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

        for t=1:reco{cpt}(3)     
            K = sigma_new*PSI_update{t}' * inv(accuracy*eye(size(PSI_update{t}*sigma_new*PSI_update{t}')) + PSI_update{t}*sigma_new*PSI_update{t}');
            mu_new = mu_new + K* (ynew{t}' - PSI_update{t}*mu_new);
            sigma_new = sigma_new - K*(PSI_update{t}*sigma_new);
  
        end

        newmu{trial}{cpt} = mu_new;
        newSigma{trial}{cpt} = nearestSPD(sigma_new);
        %[~,p] = chol(sigma_new)
        %disp(['value of p ', num2str(p)]);
        telapsed(trial,cpt) = toc(tstart);
    %     PSI_inf = computeBasisFunction (z,nbFunctions, nbDof, alpha_inf, floor(100/alpha_inf), center_gaussian, h, floor(100/alpha_inf));
    %     trajInf = (PSI_inf*newmu{cpt});
    %     minTime = min(size(trajInf,1), size(y_trial_Tot{1},1))/nbDofTot;
    %    trajErr{cpt} = abs(trajInf(1:minTime) - y_trial_Tot{1}(1:minTime));
        %reco{cpt} = [reco{cpt-1}(1),reco{cpt-1}(2),reco{cpt-1}(3)]; %{trajectory recognized, probability, timestep} of the best likelihood
    %     
    %     clear u ;
    %     u = psiTrial*mu_new;
    %     S = psiTrial*sigma_new*psiTrial' + accuracy*eye(size(psiTrial*sigma_new*psiTrial'));
    %     logkike(step) = logLikelihood( y_trial_Tot{trial}', u', S); 
        PSI_inf{trial}{cpt} = computeBasisFunction (z,nbFunctions, nbDof, alpha_inf{trial}(cpt), floor(100/alpha_inf{trial}(cpt)), center_gaussian, h, floor(100/alpha_inf{trial}(cpt)));
        cpt = cpt+1;

        clear mu_new sigma_new u psiTrial;
    end
    cpt = cpt-1;
    stock_tps(trial,:) = telapsed(trial, 1:totstep);
    trajInf{trial} = (PSI_inf{trial}{cpt}*newmu{trial}{cpt});
    minTime = min(size(trajInf{trial},1), size(y_trial_Tot{1}{trial},1))/nbDofTot;
    clear trajErr;
    trajErr = abs(trajInf{trial}(1:minTime) - y_trial_Tot{1}{trial}(1:minTime));
    stock_err{trial} = trajErr';
%end
 
%clear cpt mu_w_coord mu_w_f sigma_w_coord PSI_coor PSI_forces PSI_mean u sigma prob reco mu_n sigma_n y_trial_nbData realAlpha PSI_update ynew K mu_new sigma_new