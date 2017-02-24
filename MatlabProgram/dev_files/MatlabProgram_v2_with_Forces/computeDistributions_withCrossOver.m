% In this part, we compute the distribution for each kind of trajectories.
for i=1:nbKindOfTraj
   
    for j = 1:var(i)
        %we compute the phasis
        alpha2{i}(j) = z / totalTime(i,j);
        %we compute the corresponding PSI matrix
        PSI{i}{j} = computeBasisFunction (z,nbFunctions, nbDof, alpha2{i}(j), totalTime(i,j), center_gaussian, h, totalTime(i,j));
    end
    %alphaTest(i) = z / totalTimeTrial(i); %cross over test
    alphaTest2(i) = z / totalTimeTrial2(i); %inference test
    %PSI_test{i} = computeBasisFunction (z,nbFunctions, nbDof, alphaTest(i), totalTimeTrial(i), center_gaussian, h, totalTimeTrial(i));
    PSI_test2{i} = computeBasisFunction (z,nbFunctions, nbDof, alphaTest2(i), totalTimeTrial2(i), center_gaussian, h, totalTimeTrial2(i));
    
    mu_alpha(i) = mean([alpha2{i}, alphaTest2(i)]);
    sigma_alpha(i) = cov([alpha2{i}, alphaTest2(i)]);
    min_alpha_i(i) = min([alpha2{i}, alphaTest2(i)]);
    max_alpha_i(i) = max([alpha2{i}, alphaTest2(i)]);
end
min_alpha = min(min_alpha_i(i));
max_alpha = max(max_alpha_i(i));

PSI_z = computeBasisFunction (z,nbFunctions,nbDof, 1, z,center_gaussian,h, z);

%w computation for each trials
for i = 1 : nbKindOfTraj
        
    val = 0;
    for cpt =1:size(nbDof,2)
        val = val + nbDof(cpt)*nbFunctions(cpt);
    end
    mu_w{i} = zeros(nbTotFunctions, 1);  
    sigma_w{i} = zeros(nbTotFunctions, nbTotFunctions);  

    for nbCrossOver=1:var(i)
    tstart = tic;

        w = zeros(var(i),val);

        for j = 1 : var(i)
            %resolve a little bug
            sizeY  = size(y{i}{j},1);
            if(sizeY ~= size(PSI{i}{j},1))
                y{i}{j} = y{i}{j}(1:sizeY-(nbDofTot));
                totalTime(i,j) = totalTime(i,j) -nbDofTot;
                alpha2{i}(j) = z /totalTime(i,j);
            end
    
            w(j,:) = (PSI{i}{j}'*PSI{i}{j}+1e-12*eye(val)) \ PSI{i}{j}' * y{i}{j};        
           
        end
        %computation of the w distribution     
        muTraining{i,nbCrossOver}(:,1) = mean(w)';
    %     sigma_w{i} = zeros(size( mu_w{i},1),size( mu_w{i},1));
    %     for j=1:var(i)
    %        sigma_w{i} =sigma_w{i} + (w{i}(j,:)' -  mu_w{i})*(w{i}(j,:)' - mu_w{i})';
    %     end
    %     sigma_w{i}= sqrt(sigma_w{i} / (var(i)));
        sigmaTraining{i,nbCrossOver} = nearestSPD(cov(w)); %sometimes have < 0 for forces as it is not
        
        %linked at all
     % sigma_w{i} = nearestSPD(sigma_w{i});
        wtest = (PSI_test2{i}'*PSI_test2{i}+1e-6*eye(val)) \ PSI_test2{i}' * y_trial_Tot2{i};      
       % logLike(i,nbCrossOver) = logLikelihood(wtest,muTraining{i,nbCrossOver},  sigmaTraining{i,nbCrossOver});
        dist{i}(nbCrossOver,:) = abs(y_trial_Tot2{i}(1/alphaTest2:1/alphaTest2:totalTimeTrial2*nbDofTot) - PSI_z*muTraining{i,nbCrossOver});
%         if(nbCrossOver > 1)
%             logLike2(i,nbCrossOver-1) = logLikelihood(muTraining{i,nbCrossOver-1},muTraining{i,nbCrossOver},  sigmaTraining{i,nbCrossOver});
%         end
    
        %change the subset used to do the test
        switchData;
        mu_w{i} = mu_w{i} + muTraining{i,nbCrossOver};
        sigma_w{i} =  sigma_w{i} + sigmaTraining{i,nbCrossOver};
         tpCompute(nbCrossOver) = toc(tstart);
       %  plot(dist{i}(nbCrossOver,:));hold on;
        %clear w
    end
    mu_w{i} = mu_w{i} /var(i);
    %err(i,:) = mean(logLike(i,:));
%    for nbCrossOver=1:var(i)
%        sigma_w{i} = sigma_w{i} + (muTraining{i,nbCrossOver} - mu_w{i})* (muTraining{i,nbCrossOver} - mu_w{i})';
%    end
    sigma_w{i} = sigma_w{i} / var(i);
    %[~,p] = chol(sigma_w{i})
end

 %figure;
  %  plot(err);

   % clear i j test w PSI min_alpha_i max_alpha_i cpt nf3D sizeY type val;
