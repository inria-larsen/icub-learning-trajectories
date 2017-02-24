 clear logTest w
 %close all
 fig = figure;
subplot(3,1,1);
 PSI_inf = computeBasisFunction (z,nbFunctions, nbDof, alpha_inf, floor(100/alpha_inf), center_gaussian, h, floor(100/alpha_inf));
 PSI_real = computeBasisFunction (z,nbFunctions, nbDof, alphaTest, floor(100/alphaTest), center_gaussian, h, floor(100/alphaTest));
 PSI_mean = computeBasisFunction (z,nbFunctions, nbDof, mu_alpha, floor(100/mu_alpha), center_gaussian, h, floor(100/mu_alpha));
 % PSI_inf  = PSI_mean;
 w = (PSI_real'*PSI_real+1e-6*eye(size(PSI_real'*PSI_real))) \ PSI_real' * y_trial_Tot{trial}(1:nbDofTot*floor(100/alphaTest)); 
 
 for(i=1:size(newmu,2))
    logTest(i) = logLikelihood(w, newmu{i}, newSigma{i});
    fig(size(fig,2) + 1) = plot((i-1)*floor(totalTimeTrial(trial) /10)+1, logTest(i), 'x'); hold on;
 end

fig(size(fig,2) + 1) = plot([1:floor(totalTimeTrial(trial) /10):floor(totalTimeTrial(trial) /10)*size(newmu,2)], logTest, '--'); hold on;
a = size(fig,2);
title('Loglikelihood according to the actual update');
xlabel('Iterations');
ylabel('log(p)');

subplot(3,1,2);

courb = PSI_mean*mu_w{1};
fig(size(fig,2) + 1) = plot(courb(1:z/mu_alpha), 'r');hold on
b =size(fig,2);

%for (i=2:size(newmu,2))
i=size(newmu,2)

courb = PSI_inf*newmu{i};
fig(size(fig,2) + 1) = plot(courb(1:z/alpha_inf), 'b');hold on;
%fig(size(fig,2) + 1) = plot(PSI_real*newmu{i}, 'm');hold on;
%end
fig(size(fig,2) + 1) = plot(y_trial_Tot{trial}(1:totalTimeTrial(trial)), 'g');hold on;



title('Infered trajectory with velocity infered');
xlabel('Iterations');
ylabel('input value');
legend('distribution learnt', 'distribution updated', 'real trial');

subplot(3,1,3);


trajInf = (PSI_inf*newmu{i});

minTime = min(size(trajInf(1:totalTimeTrial(trial)),1), size(y_trial_Tot{trial},1))/nbDofTot;
plot(abs(trajInf(1:minTime) - y_trial_Tot{trial}(1:minTime)));
    annotation('textbox', [.2 .15 .1 .1], 'String', ['mean(error) =', num2str(mean(abs(trajInf(1:minTime) - y_trial_Tot{trial}(1:minTime))))]); 

title('Error representation of the infered trajectory')
xlabel('iterations')
ylabel('abs(y_{real} - y_{infered})')
