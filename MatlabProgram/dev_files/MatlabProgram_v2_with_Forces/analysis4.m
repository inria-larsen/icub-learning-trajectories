%Here we draw all the plot for the meeting
i=1 ;% type de trajectoires
    set(0,'DefaultAxesFontSize',18)
close all
fig = figure;

minTime= min([ totalTime ,totalTimeTrial]);
subplot(2,2,3); %err representatoon according to time
for i=1:10
    for cpt=1:size([floor( minTime/nbStep): floor(minTime/nbStep):  minTime],2)
        
        trajInf{i}{cpt} = PSI_inf{i}{cpt}*newmu{i}{cpt};
        trajIntTime = size(trajInf{1}{1},1);
        err(i,cpt) = 0;
        for dof=1:nbDofTot
            err(i,cpt) = err(i,cpt) + abs(trajInf{i}{cpt}(trajIntTime - dof+1) - y_trial_Tot{1}{i}(totalTimeTrial(i) - dof +1));
        end
    end
end
    boxplot(err, [floor( minTime/nbStep): floor(minTime/nbStep):  minTime] );
    title('Final position error')
    xlabel('# Samples known');
    ylabel('\Sigma_d_i_m (errors)') ;

subplot(2,2,4);
boxplot(stock_tps, [floor( minTime/nbStep): floor(minTime/nbStep):  minTime] );
title('Time to infer');
xlabel('# samples known')
%yaxis(0.0395, 0.0415)

i=1;
exp = 2
subplot(2,2,2);
visualisationShared(PSI_z*mu_w{i}, PSI_z*1.96*sqrt(diag(sigma_w{i} )), nbDofTot, z,  1, 'b', fig);hold on;
fig(size(fig,2)+1) = plot(PSI_z*mu_w{i});
distr = size(fig,2);
fig(size(fig,2)+1) = plot(y_trial_Tot{i}{exp}(1:totalTimeTrial(exp)/z:totalTimeTrial(exp)), 'r');
dataa = size(fig,2);
fig(size(fig,2)+1) = plot(y_trial_Tot{i}{exp}(1:totalTimeTrial(exp)/z:z/floor(minTime/nbStep)*3*nbDofTot), '+r');
for  cpt=1:size([floor( minTime/nbStep): floor(minTime/nbStep):  minTime],2)
    if(cpt==3)
        fig(size(fig,2)+1) = plot(PSI_z*newmu{exp}{cpt}, 'm');
        threee = size(fig,2);
    else
        fig(size(fig,2)+1) = plot(PSI_z*newmu{exp}{cpt}, ':', 'Color', [0, cpt/size([floor( minTime/nbStep): floor(minTime/nbStep):  minTime],2),0 ]);
    end
end
title('Example of an infered trajectory');


legend(fig([distr dataa (dataa+1) threee (threee+1)]), 'Distribution', 'data', 'data known', 'inference 24 samples', 'other inferences')