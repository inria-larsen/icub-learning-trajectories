 clear logTest w
 %close all
trial = 1
    set(0,'DefaultAxesFontSize',18)
 fig = figure;
subplot(2,2,1);


% PSI_inf = computeBasisFunction (z,nbFunctions, nbDof, alpha_inf, floor(100/alpha_inf), center_gaussian, h, floor(100/alpha_inf));
 PSI_real = computeBasisFunction (z,nbFunctions, nbDof, 100/ totalTimeTrial(1,trial), totalTimeTrial(1,trial), center_gaussian, h, totalTimeTrial(1,trial));
 PSI_mean = computeBasisFunction (z,nbFunctions, nbDof, mu_alpha, floor(100/mu_alpha), center_gaussian, h, floor(100/mu_alpha));
  PSI_inf  = PSI_mean;
  alpha_inf = mu_alpha
 w = (PSI_real'*PSI_real+1e-6*eye(size(PSI_real'*PSI_real))) \ PSI_real' * y_trial_Tot{1}{trial}(1:nbDofTot*totalTimeTrial(1,trial)); 
 
 for(i=1:size(newmu{trial},2))
    logTest(i) = logLikelihood(w, newmu{trial}{i}, newSigma{trial}{i});
    fig(size(fig,2) + 1) = plot((i-1)*floor(totalTimeTrial(1,trial) /10)+1, logTest(i), 'x'); hold on;
 end

fig(size(fig,2) + 1) = plot([1:floor(totalTimeTrial(1,trial) /10):floor(totalTimeTrial(1,trial) /10)*size(newmu{trial},2)], logTest, '--'); hold on;
a = size(fig,2);
title('Loglikelihood according to the actual update');
xlabel('Iterations');
ylabel('log(p)');

subplot(2,2,2);

courb = PSI_mean*mu_w{1};
fig(size(fig,2) + 1) = plot(courb(1:z/mu_alpha), 'r');hold on
b =size(fig,2);

%for (i=2:size(newmu{trial},2))
i=size(newmu{trial},2)

courb = PSI_inf*newmu{trial}{i};
fig(size(fig,2) + 1) = plot(courb(1:z/alpha_inf), 'b');hold on;
%fig(size(fig,2) + 1) = plot(PSI_real*newmu{trial}{i}, 'm');hold on;
%end
fig(size(fig,2) + 1) = plot(y_trial_Tot{1}{trial}(1:totalTimeTrial(1,trial)), 'g');hold on;



title('Infered trajectory with velocity infered');
xlabel('Iterations');
ylabel('input value');
legend('distribution learnt', 'distribution updated', 'real trial');
annotation('textbox', [.25 .7 .1 .1], 'String', ['alpha_{expected} =', num2str(alpha_inf), ' alpha_{real} =', num2str(mu_alpha)], 'FontSize', 16); 

subplot(2,2,3);


trajInf = (PSI_inf*newmu{trial}{i});

minTime = min(size(trajInf,1), size(y_trial_Tot{1}{trial},1))/nbDofTot;
plot(abs(trajInf(1:minTime) - y_trial_Tot{1}{trial}(1:minTime)));
te =  annotation('textbox', [.2 .1 .1 .1], 'String', ['mean(error) =', num2str(mean(abs(trajInf(1:minTime) - y_trial_Tot{1}{trial}(1:minTime))))], 'FontSize', 16); 

te.FontSize = 18;
title('Error representation of the infered trajectory')
xlabel('iterations')
ylabel('abs(y_{real} - y_{infered})')


subplot(2,2,4)
plot(telapsed, '*');
title('Time elpased according to the number of iterations used to do the inference');
xlabel('Iterations used to do the inference');
ylabel('Coputation time (s)');


 %%%%%%%%%%%%%%%%INFERENCE ANALYSIS ACCORDING TO TIME
close all
fig3 = figure(3);
filename = 'inference.gif';
i=1; type =1;
trial=1;
vaCl=PSI_inf{trial}{size(newmu{trial},2)}*mu_w{i}; 
plot(vaCl((1:totalTimeTrial(1,trial))),'r');hold on;
vacl2 = PSI_inf{trial}{size(newmu{trial},2)}*1.96*sqrt(diag(sigma_w{i}));
visualisationShared(PSI_z*mu_w{i}, PSI_z*1.96*sqrt(diag(sigma_w{i})), nbDofTot, z,type, 'k', fig3);
plot(y_trial_Tot{1}{trial}(1:totalTimeTrial(trial)), '.m');
 for value=1:var(i)
    plot(y{i}{value}(1: floor(totalTime(1,value) / z):z), ':k');hold on;
 end
    drawnow
    frame = getframe(3);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
for j = 1:size(newmu{trial},2)-1
    vaCl =PSI_inf{trial}{j}*newmu{trial}{j};
    h = plot(vaCl(1:totalTimeTrial(1,trial)), 'Color', [0,j/size(newmu{trial},2),0]);
    %hh= visualisationShared(PSI_inf*newmu{j}, PSI_z*1.96*sqrt(diag(newSigma{j})), nbDofTot, z,type, 'g', figAnalysis);
    h2 = plot(floor(totalTimeTrial(trial) /10)*j, y_trial_Tot{i}{trial}(floor(totalTimeTrial(trial) /10)*j), '*', 'Color', [0,j/size(newmu{trial},2),0]);
    drawnow
    frame = getframe(3);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,filename,'gif','WriteMode','append');
    delete(h); 
    %delete(hh);
end
vaCl =PSI_inf{trial}{size(newmu{trial},2)}*newmu{trial}{size(newmu{trial},2)}; 
vacl2 = PSI_inf{trial}{size(newmu{trial},2)}*1.96*sqrt(diag(newSigma{trial}{size(newmu{trial},2)}));
 plot(vaCl(1:totalTimeTrial(1,trial)), 'g');
 visualisationShared(vaCl, vacl2, nbDofTot, z,type, 'g', fig3);
 xlabel('Iterations');
    ylabel('data value');
    title('inference according to time');
     drawnow
    frame = getframe(3);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,filename,'gif','WriteMode','append');
    
   
%     
