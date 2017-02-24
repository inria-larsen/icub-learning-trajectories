%in this function, we try to recongize a movement from the 30 first data
%and to complete it.
%TO do that, we consider the phasis of the movement as the mean of the
%phasis used during the learning.

  %TODO A SUPPR pour correspondre aux data tests
    y_trial_Tot = y_trial_tot;
   %totalTimeTrial = totalTime;
    nbDofTot = 1;
   
    %accuracy that we suppose the data has 
    accuracy=0.0001;
    
for trial=1:10
%trial=2

    %computation of the loglikelihood for each trajectory using only cartesian
    %coordinates
    typeReco = 1; % what kind of input we use for the reco{cpt}
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

    %%TODO will be suppr (to draw reco according to time)
    nbStep= 5;
    cpt = 1;
     fig = figure  
        
    for valstep=floor( min([ totalTime ,totalTimeTrial]) /nbStep): floor(min([ totalTime ,totalTimeTrial])/nbStep):  min([ totalTime ,totalTimeTrial])
        
%4     
        tstart = tic;
        reco{cpt} = [0.0 , 0.0, 0.0 ];
        
        for i=1:nbKindOfTraj
            %we compute the probability it correspond to the actual trial of
            %"step" steps.
            y_trial_part{valstep} = [];
            for k=1:nbDof(typeReco)
               y_trial_part{valstep} = [y_trial_part{valstep} ; y_trial_Tot{1}{trial}(totalTimeTrial(1, trial)*(prevDof+k -1) + 1 : totalTimeTrial(1,trial)*(prevDof+k -1) + valstep)];
            end
            
            %%%TODO will be deleted to verify perf
            alphaTest = 100/totalTimeTrial(trial);
            inf2(1) = 1; 
            inf2(2)= alpha2{i}(1);
            %%%END TODO
            %figure(valstep)
            for numAlpha =1:var(i) 
       
                %compute the mask to know the data we known;
                 ma =  ones(1,valstep);
                 mb = zeros(1,totalTime(i,numAlpha)-valstep); 
                 ma = [ma, mb]; % 1 for data we know 0 for the others
                 mk = []; % C/C the mask for other data
                 for vv = 1:nbDofTot
                     mk = [mk, ma]; 
                 end 
                 mask{numAlpha} = logical(mk);
                 %%Choose velocity according to air under the courb
                 courb{numAlpha} = abs(PSI{i}{numAlpha}(mask{numAlpha},:)*mu_w_reco{i} - y_trial_part{valstep});                
                 distAlpha(numAlpha,:) = sum(courb{numAlpha});
                %%Not efficient
            % compute w : (PSI{i}{j}'*PSI{i}{j}+1e-12*eye(val)) \ PSI{i}{j}' * y{i}{j}
           	PSI_reco_alpha{valstep}{i} = computeBasisFunction(z,nbFunctions(typeReco), nbDof(typeReco), alpha2{i}(numAlpha), z/alpha2{i}(numAlpha), center_gaussian(typeReco), h(typeReco),valstep);

            wtest = (PSI_reco_alpha{valstep}{i}'*PSI_reco_alpha{valstep}{i}+ 1e-12*eye(size(PSI_reco_alpha{valstep}{i}'*PSI_reco_alpha{valstep}{i},1))) \ PSI_reco_alpha{valstep}{i}'* y_trial_part{valstep};    
            %proba traj & alpha
            probA(numAlpha) = logLikelihoodWithAlpha(wtest,mu_w_reco{i},sigma_w{i}, alpha2{i}(numAlpha), mu_alpha, sigma_alpha);
           % disp(['proba alpha= ', num2str(logLikelihood(alpha2{i}(numAlpha), mu_alpha, sigma_alpha)), ' proba traj =', num2str(probA(numAlpha))]);
            probB(numAlpha) = logLikelihood(wtest,mu_w_reco{i},sigma_w{i});
            errorAlpha(numAlpha) = abs(alpha2{i}(numAlpha) - (100/totalTimeTrial(trial)));
%             %%%TODO will be deleted to verify perf
%                 if(abs(inf2(2) - alphaTest) > abs(alpha2{i}(numAlpha) -  alphaTest))
%                     inf2(1) = numAlpha;
%                     inf2(2) =  alpha2{i}(numAlpha);
%                 end
%                 %%%END TODO 
                
%%%ANOTHER TEST            
%             subplot(4,1,3);
%             scatter(numAlpha,probA(numAlpha));hold on;%proba traj
%             title('proba traj et alpha');
%             subplot(4,1,1)
%             scatter(numAlpha, logLikelihood(wtest,mu_w_reco{i},sigma_w{i}));hold on;
%             title('proba traj');
%             subplot(4,1,2)
%             scatter(numAlpha, logLikelihood(alpha2{i}(numAlpha), mu_alpha, sigma_alpha));hold on;
%             title('proba alpha');
%             subplot(4,1,4)
%             scatter(numAlpha, abs(alpha2{i}(numAlpha) - (100/totalTimeTrial(trial))));hold on;
%             title('real error alpha');
%%%END TEST

            end
            
            [valMin, indMin] = min(distAlpha);
            
            %%TODO CHEAT TEST ALPHA
            indMin = trial;
            
            
            
%             [valMin, indMin2] = max(probA);%min(distAlpha);
%             [valErr, indErr] = min(errorAlpha);
%             if(indMin - indErr ~=0)
%                 disp(['err for the ' num2str(valstep), ' we have ', num2str(indMin), ' should be ',num2str(indErr) ]);
%                 figure;
%                 plot(y_trial_part{valstep}, 'b');hold on;
%                 plot(PSI{i}{indMin}(mask{indMin},:)*mu_w_reco{i}, 'm');
%                 plot(PSI{i}{indMin2}(mask{indMin2},:)*mu_w_reco{i}, 'k');
%                 plot(PSI{i}{indErr}(mask{indErr},:)*mu_w_reco{i}, 'g');
%                 legend('real traj', 'maxLogWithAlpha', 'maxLog', 'minErr');
%                 
%             else
%                 disp(['inf traj ok ind=', num2str(indErr)]);
%             end
            

                %%%TODO CHHHHEEEAAT
             alpha_inf{trial}(cpt) = 100/totalTimeTrial(trial)%alpha2{i}(indMin);
            
             %TODO suppr these two lines
             %[~,indalpha] = min(abs(alpha2{i}(numAlpha) - (100/totalTimeTrial(trial)))); % cheat to verify the computation
            % disp(['We supposed ', num2str(alpha_inf{trial}(cpt)), ' instead of ', num2str(inf2(2)),' for the real: ', num2str(alphaTest)]);
                    
%                     figure;
%                      subplot(2,2,1);
%                      plot(y_trial_part{valstep}, 'b');hold on;
%                      plot(PSI{i}{indMin}(mask{indMin},:)*mu_w_reco{i}, 'm');
%                      plot(PSI{i}{inf2(1)}(mask{inf2(1)},:)*mu_w_reco{i}, 'g');
%                      legend('reel', 'quon a inferé', 'avec le alpha reelment plus proche');
%                      subplot(2,2,2);
%                      plot(courb{indMin},'m');hold on;
%                      plot(courb{inf2(1)},'g');
%                      legend('we have', 'we should have less');  
             
            PSI_reco{valstep}{i} = computeBasisFunction(z,nbFunctions(typeReco), nbDof(typeReco), alpha_inf{trial}(cpt), floor(z/alpha_inf{trial}(cpt)), center_gaussian(typeReco), h(typeReco),valstep);
            PSI_tot{valstep}{i} = computeBasisFunction(z,nbFunctions, nbDof, alpha_inf{trial}(cpt), floor(z/alpha_inf{trial}(cpt)), center_gaussian, h,valstep);
            
            %%Not efficient
            % compute w : (PSI{i}{j}'*PSI{i}{j}+1e-12*eye(val)) \ PSI{i}{j}' * y{i}{j}
            wtest = (PSI_reco{valstep}{i}'*PSI_reco{valstep}{i}+ 1e-12*eye(size(PSI_reco{valstep}{i}'*PSI_reco{valstep}{i},1))) \ PSI_reco{valstep}{i}'* y_trial_part{valstep};    
            prob(i,cpt) = logLikelihoodWithAlpha(wtest,mu_w_reco{i},sigma_w_reco{i}, alpha_inf{trial}(cpt), mu_alpha, sigma_alpha);
           
            %%TODO prendre ça plutot
            %%we take the air between courbs whearas the loglikelihood
            %dist(i,cpt) = mean(abs(y_trial_part{valstep} - PSI_reco{valstep}{i}*mu_w_reco{i}));
            
            %we record the max of probability to know wich distribution we
            %recognize
            reco{cpt}(2) = prob(i,cpt);%dist(i,cpt)
            reco{cpt}(1) = i;
            reco{cpt}(3) = valstep;
        end

        %we retrieve the computed distribution that correspond to the recognized
        %trajectory
        mu_new = mu_w{reco{cpt}(1)};
        sigma_new = sigma_w{reco{cpt}(1)}; 

        %%we complete the data with others dimension of the movement
        y_prev = [];
        %first, dim before the one(s) used for the recognition
        for i=1:prevDof
            y_prev  = [y_prev ; y_trial_Tot{1}{trial}(totalTimeTrial(1, trial)*(i -1) + 1 : totalTimeTrial(1,trial)*(i -1) + reco{cpt}(3))];
        end
        %then, dim after the one(s) used for the reco
        y_after = [];
        for i= (1 + prevDof + nbDof(typeReco)): nbDofTot
           y_after = [y_after;  y_trial_Tot{1}{trial}(totalTimeTrial(1, trial)*(i -1) + 1 : totalTimeTrial(1, trial)*(i -1) + reco{cpt}(3))];   
        end
        %all data (all DOF) for time t=1:nbData
        y_trial_nbData =[y_prev ; y_trial_part{reco{cpt}(3)}; y_after ];

%%%%STAT        
        %we aren't suppose to know "realData",  it is only used to draw the real
        %trajectory of the sample if we continue it to the end
        realAlpha = z /totalTimeTrial(1, reco{cpt}(1));
%%%%END_STAT
        
%%%%%%%%%%%%%%OLD VERSION        
%we need to have the psi matrix and vector value according to time to
%update the distribution (just a rewriting of data to simplify the next
%computation.
%         vr=0;
%         for j=1:typeReco-1
%             vr = vr + nbDof(j);
%         end
%         for t=1:reco{cpt}(3) %reco{cpt}(3) = valstep number of step we have
%             for i=vr + 1: vr + nbDof(typeReco)
%                 PSI_update{t}(i,:) = PSI_tot{valstep}{reco{cpt}(1)}(t + reco{cpt}(3)*(i-1),:);
%                 ynew{t}(i) = y_trial_nbData(t + reco{cpt}(3)*(i-1)) ;
%             end
%         end
% 
%         for t=1:reco{inf.figcpt}(3) %=valStep    
%             K = sigma_new*PSI_update{t}' * inv(accuracy*eye(size(PSI_update{t}*sigma_new*PSI_update{t}')) + PSI_update{t}*sigma_new*PSI_update{t}');
%             mu_new = mu_new + K* (ynew{t}' - PSI_update{t}*mu_new);
%             sigma_new = sigma_new - K*(PSI_update{t}*sigma_new);
%         end
%%%%%%%%%%%%%%END_OLD_VERSION        


        %%Creation of the basis function with good velocity & nbData
        %creation of the mask to have only basis linked to the data known
        ma =  ones(1,valstep);
        mb = zeros(1,totalTime(i,indMin)-valstep); 
        ma = [ma, mb]; % 1 for data we know 0 for the others
        mk = []; % C/C the mask for other data
        for vv = 1:nbDofTot
        	mk = [mk, ma]; 
        end 
        mask{indMin} = logical(mk);
        %creation  of the basis function matrix
        PSI_update = PSI{reco{cpt}(1)}{indMin}(mask{indMin},:);

       %update
       K = sigma_new*PSI_update' * inv(accuracy*eye(size(PSI_update*sigma_new*PSI_update')) + PSI_update*sigma_new*PSI_update');
       mu_new = mu_new + K* (y_trial_nbData - PSI_update*mu_new);
       sigma_new = sigma_new - K*(PSI_update*sigma_new);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%STAT   
subplot(3,2,cpt);
        newmu{trial}{cpt} = mu_new;
        newSigma{trial}{cpt} = nearestSPD(sigma_new);
        telapsed(trial,cpt) = toc(tstart);
        PSI_inf{trial}{cpt} = PSI{reco{cpt}(1)}{indMin};
        % = computeBasisFunction (z,nbFunctions, nbDof, alpha_inf{trial}(cpt), floor(100/alpha_inf{trial}(cpt)), center_gaussian, h, floor(100/alpha_inf{trial}(cpt)));
        %%%PLOT  
        plot(y_trial_nbData, '+b');hold on;
        plot(y_trial_Tot{1}{trial}, ':b');
        plot(PSI_inf{trial}{cpt}*newmu{trial}{cpt}, 'Color', [cpt/10,0,0]);
        ba = visualisationShared(PSI_inf{trial}{cpt}*mu_w{reco{cpt}(1)}, PSI_inf{trial}{cpt}*1.96*sqrt(diag(sigma_w{reco{cpt}(1)})), nbDofTot, totalTime(reco{cpt}(1),indMin),1, 'g', fig);
        plot(PSI_inf{trial}{cpt}*mu_w{1}, 'g');
        name =['Infered trajectory with ', num2str(valstep), ' iterations known'];
        if(cpt==3)
            newax = subplot(3,2,6);
            ba = visualisationShared(PSI_z*mu_w{reco{cpt}(1)}, PSI_z*1.96*sqrt(diag(sigma_w{reco{cpt}(1)})), nbDofTot, z,1, 'm', newax);hold on;
            a = size(newax,2);
            newax(size(newax,2)) = plot(PSI_z*mu_w{1}, 'm');  
            ba = visualisationShared(PSI_z*mu_new, PSI_z*1.96*sqrt(diag(sigma_new)), nbDofTot, z,1, 'g', newax);
            newax(size(newax,2) + 1) = plot(PSI_z*mu_new, 'g');  
            legend(newax([a size(newax,2)]), 'previous distribution' , 'new distribution');
            title(['distribution before after the update with ', num2str(valstep), 'samples']);
            xlabel('Iterations');
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END_STAT
        cpt = cpt+1;
    end
    
%%%%PLOT    
%         title(name);
%         xlabel('Iterations');
%         ylabel('trajectorie data');
%         legend('samples known' , 'total Trajectory', 'infered trajectory', 'initial learned distribution');

%%STAT
    cpt = cpt-1;
    stock_tps(trial,:) = telapsed(trial, 1:5);
    %trajInf{trial} = (PSI_inf{trial}{cpt}*newmu{trial}{cpt});
    %minTime = min(size(trajInf{trial},1), size(y_trial_Tot{1}{trial},1))/nbDofTot;
    %clear trajErr;
    %trajErr = abs(trajInf{trial}{step}(1:minTime) - y_trial_Tot{1}{trial}(1:minTime));
    %stock_err{trial} = trajErr';
end