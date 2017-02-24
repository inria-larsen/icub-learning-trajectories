%Here we draw all the plot for the meeting
i=1 ;% type de trajectoires

close all;

% 
% %%%%%In figure one all information about the computation of the
% %%%%%distribution

fig = figure(1);
subplot(2,2,2); 
fig(size(fig,2)+1) = visualisationShared(PSI_z*mu_w{i}, PSI_z*1.96*sqrt(diag(sigma_w{i} )), nbDofTot, z,  1, 'b', fig)
for value=1:var(i)
    subplot(2,2,1);
    fig(size(fig,2)+1) = plot(y{i}{value}, 'c');hold on;
    subplot(2,2,2);    
   % numPlot = PSI_z*muTraining{i,value};
    fig(size(fig,2)+1) = plot(numPlot(1:100), 'Color', [1, 0, 0]); hold on;
    title('Distribution learned from the cross validation processing');
    %p(value) = legend(['crossover n°' num2str(value)])
    subplot(2,2,3); 
    fig(size(fig,2)+1) = plot(dist{i}(value,1:100), 'Color', [1, 0, 0]); hold on
end
fig(size(fig,2)+1) = plot(mean(dist{i}(:,1:100)), 'b');
b = size(fig,2);
title('Test during the cross validation (difference between courbs)');
xlabel('Iterations');
ylabel('error data');
subplot(2,2,1);
xlabel('Iterations');
ylabel('data value');
title('data samples for learning (blue) and the trial (green)');

fig(size(fig,2)+1) = plot(y_trial_Tot{1}{i}(1:100/totalTimeTrial(1):totalTimeTrial(1)), 'g');hold on;
a =size(fig,2);
%fig(size(fig,2)+1) = plot(y_trial_Tot2{1}{i}(1:alphaTest2:totalTimeTrial2(1)), 'c');hold on;
%legend(fig([size(fig,2)]), 'trial we keep to test the inference');

subplot(2,2,2);
numPlot = PSI_z*mu_w{i};
fig(size(fig,2)+1) = plot(numPlot(1:100), 'b');
c =size(fig,2);
xlabel('Iterations');
ylabel('data value');
%xaxis([70 100]);

subplot(2,2,4)
fig(size(fig,2) + 1) = plot(tpCompute, '+r');hold on;
fig(size(fig,2) + 1)= plot([1:1:size(tpCompute,2)],mean(tpCompute)*ones(1,size(tpCompute,2)));hold on;
d = size(fig,2);
xlabel('cross validation step');
ylabel('computation time');
title('Time used to compute each step of the cross validation');

legend(fig([a b c d]), 'data trial', 'mean error', 'mean distribution', 'mean computation time')

%%%%%%%%%%inference analysis %%%%%%%%%%%%%%%%%
% 
fig = figure(2);
%nbTrial = size(stock_err,2);
clear courb difMu
%%%%%%%%%%loglikelihood
subplot(2,2,1);
PSI_mean = computeBasisFunction (z,nbFunctions, nbDof, mu_alpha, floor(100/mu_alpha), center_gaussian, h, floor(100/mu_alpha));
difmean = zeros(nbStep);
for trial=1:nbTrial
    %PSI_inf = computeBasisFunction (z,nbFunctions, nbDof, alpha_inf, floor(100/alpha_inf), center_gaussian, h, floor(100/alpha_inf));
    PSI_real{trial} = computeBasisFunction (z,nbFunctions, nbDof, 100 / totalTimeTrial(trial), totalTimeTrial(trial), center_gaussian, h, totalTimeTrial(trial));

    %PSI_inf  = PSI_mean;
    %alpha_inf = mu_alpha
    wtrial{trial} = (PSI_real{trial}'*PSI_real{trial}+1e-12*eye(size(PSI_real{trial}'*PSI_real{trial}))) \ PSI_real{trial}' * y_trial_Tot{1}{trial}(1:nbDofTot*totalTimeTrial(trial));
   % figure(13+trial);
   % plot(wtrial{trial},'r');hold on;
    
    for(i=1:size(newmu{trial},2)) % for each step (more known data at each step)
    %    plot(abs(newmu{trial}{i}- wtrial{trial}), '-+', 'Color', [0, 0, i/size(newmu{trial},2)]);hold on;
        difMu(trial,i) = mean(abs(wtrial{trial}- newmu{trial}{i}));%logLikelihood(wtrial{trial}, newmu{trial}{i}, sigma_w{1});%newSigma{trial}{i})
    end
    %difmean = difmean + difMu(trial,:);
 %   figure(2)
  %  subplot(2,2,1);
    fig(size(fig,2) + 1) = plot([floor(totalTimeTrial(1,trial) /nbStep): floor(totalTimeTrial(1,trial) /nbStep): z/ max_alpha], difMu(trial,1:size(newmu{trial},2)), '+', 'Color', [0, 0, trial/nbTrial]); hold on;    

    fig(size(fig,2) + 1) = plot([floor(totalTimeTrial(1,trial) /nbStep): floor(totalTimeTrial(1,trial) /nbStep): z/ max_alpha], difMu(trial, 1:size(newmu{trial},2)), ':', 'Color', [0, 0, trial/nbTrial]); hold on;    
end
fig(size(fig,2) + 1) = plot( [ (floor(mean(totalTimeTrial) /4)) : (floor(mean(totalTimeTrial) /4)): floor(mean(totalTimeTrial))], mean(difMu(:,1:4))); hold on;    


 %plot([1:ceil(totalTimeTrial(trial) /size(difMu(trial,:),2)):totalTimeTrial(trial)],mean(difMu));
i=1

title('Parameters error according to the actual update');
xlabel('Iterations');
ylabel('|w_d - w_i_n_f|');
%xaxis([11 100]);


%%%%%%%%Infered trajectory
trial = 1;
subplot(2,2,2);
type =1;
courb = PSI_inf{trial}{nbStep-1}*mu_w{i};
fig(size(fig,2) + 1) = plot(courb, 'r');hold on;
%visualisationShared(courb, PSI_inf{trial}{10}*1.96*sqrt(diag(sigma_w{i})), nbDofTot, 100/ alpha_inf{trial}(10),type, 'r', type);
a =size(fig,2);

%  for value=1:var(i)
%      fig(size(fig,2)+1) = plot(y{i}{value}(1: alpha2{i}(value):totalTime(i,value)), ':k');hold on;
%  end
% d =size(fig,2);
for (i=2:size(newmu{trial},2)-1)
    %i=size(newmu,2)
  
    courb = PSI_inf{trial}{i}*newmu{trial}{i};
    fig(size(fig,2) + 1) = plot(courb, 'Color', [0, 0, i/size(newmu{trial},2)]);hold on;
    %fig(size(fig,2) + 1) = plot(PSI_real*newmu{i}, 'm');hold on;
    fig(size(fig,2) + 1) = scatter(i*floor(totalTimeTrial(1,trial) /10)+1, y_trial_Tot{1}{trial}(i*floor(totalTimeTrial(1,trial) /10)+1), '+k');hold on;
end
courb = PSI_inf{trial}{size(newmu{trial},2)}*newmu{trial}{size(newmu{trial},2)};
fig(size(fig,2) + 1) = plot(courb, '.-', 'Color', [0, 0, 1]);hold on;

c = size(fig,2);
fig(size(fig,2) + 1) = plot(y_trial_Tot{1}{trial}, '.-g');hold on;
b = size(fig,2);


title('Infered trajectory with velocity infered');
xlabel('Iterations');
ylabel('input value');

%xaxis([50 100])
legend(fig([a (a+1) c b]), 'distribution learnt', 'first update', 'last update', 'real trial');
annotation('textbox', [.25 .7 .1 .1], 'String', ['alpha_{expected} =', num2str(alpha_inf{trial}(cpt)), ' alpha_{real} =', num2str(100/totalTimeTrial(trial))], 'FontSize', 16);




%%%%%%%%%%%%%%%%%error (integral under the courbà
subplot(2,2,3);


%trajInf = (PSI_inf{trial}*newmu{trial}{i});
tmin = size(stock_err{1},2);
stockTot = stock_err{1}(1:tmin);
nbTrial = size(stock_err,2);
for i=2:nbTrial
    if(size(stock_err{i},2)< tmin)
        tmin = size(stock_err{i},2);
    end
    plot(stock_err{i}, ':');hold on;
    stockTot(1:tmin) = stockTot(1:tmin) + stock_err{i}(1:tmin);
end
stockTot = stockTot / nbTrial;
plot(stockTot(1:tmin));
minTime = min(size(trajInf,1), size(y_trial_Tot{1}{trial},1))/nbDofTot;

%plot(mean(stock_err))%abs(trajInf(1:minTime) - y_trial_Tot{trial}(1:minTime)));
annotation('textbox', [.2 .1 .1 .1], 'String', ['mean(error) =', num2str(mean(abs(trajInf{trial}(1:minTime) - y_trial_Tot{1}{trial}(1:minTime))))], 'FontSize', 16);

title('Error representation of the infered trajectory')
xlabel('iterations')
ylabel('abs(y_{real} - y_{infered})')
%yaxis([0 0.005]);


subplot(2,2,4);
i=1;
%mean(stock_tps);
%fig(size(fig,2) + 1) =plot([floor(totalTimeTrial(trial) /10):floor(totalTimeTrial(trial) /10):totalTimeTrial(trial)], mean(stock_tps), '.-r');hold on;
for j=1:nbStep
    fig(size(fig,2) + 1) =plot([floor(totalTimeTrial(trial) /nbStep):floor(totalTimeTrial(trial) /nbStep):( z/ max_alpha) ], stock_tps(j,:), ':', 'Color', [0,0,j/nbStep]);hold on;
end
title('time elapsed to infer the movement');
xlabel('number of iterations known to do the inference');
ylabel('computation time (s)');
%yaxis([2 8]);
% 
% % %%%%%%%%%%%%%%%%%INFERENCE ANALYSIS ACCORDING TO TIME
% close all
% fig3 = figure(3);
% filename = 'inference.gif';
% i=1; type =1;
% trial=4;
% vaCl=PSI_inf{trial}{size(newmu{trial},2)}*mu_w{i}; 
% plot(vaCl,'r');hold on;
% vacl2 = PSI_inf{trial}{size(newmu{trial},2)}*1.96*sqrt(diag(sigma_w{i}));
% visualisationShared(PSI_z*mu_w{i}, PSI_z*1.96*sqrt(diag(sigma_w{i})), nbDofTot, z,type, 'k', fig3);
% plot(y_trial_Tot{1}{trial}(1:totalTimeTrial(trial)), '.m');
%  for value=1:var(i)
%     plot(y{i}{value}(1: floor(totalTime(1,value) / z):z), ':k');hold on;
%  end
%     drawnow
%     frame = getframe(3);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
% for j = 1:size(newmu{trial},2)-1
%     vaCl =PSI_inf{trial}{j}*newmu{trial}{j};
%     h = plot(vaCl, 'Color', [0,j/size(newmu{trial},2),0]);
%     %hh= visualisationShared(PSI_inf*newmu{j}, PSI_z*1.96*sqrt(diag(newSigma{j})), nbDofTot, z,type, 'g', figAnalysis);
%     h2 = plot(floor(totalTimeTrial(trial) /10)*j, y_trial_Tot{i}{trial}(floor(totalTimeTrial(trial) /10)*j), '*', 'Color', [0,j/size(newmu{trial},2),0]);
%     drawnow
%     frame = getframe(3);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     imwrite(imind,cm,filename,'gif','WriteMode','append');
%     delete(h); 
%     %delete(hh);
% end
% vaCl =PSI_inf{trial}{size(newmu{trial},2)}*newmu{trial}{size(newmu{trial},2)}; 
% vacl2 = PSI_inf{trial}{size(newmu{trial},2)}*1.96*sqrt(diag(newSigma{trial}{size(newmu{trial},2)}));
%  plot(vaCl, 'g');
%  visualisationShared(vaCl, vacl2, nbDofTot, z,type, 'g', fig3);
%  xlabel('Iterations');
%     ylabel('data value');
%     title('inference according to time');
%      drawnow
%     frame = getframe(3);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     imwrite(imind,cm,filename,'gif','WriteMode','append');
%     
%    
% %     


