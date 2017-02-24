%Plot the total trial and the data we have
%nameFig = figure;

%if you want to plot on the same plot than the learnt distribution.
nameFig = figure(trial + nbDof(1));

for vff=1:nbDofTot
    subplot(nbDof(1),size(nbDof,2),vff);
    nameFig = visualisation2(y_trial_Tot{trial},sum(nbDof), totalTimeTrial(trial),reco{1}, ':m', realAlpha, nameFig);hold on;
    dtG = size(nameFig,2);
    nameFig(size(nameFig,2) + 1) = plot(y_trial{trial}(1:timeInf/z:nbData),'om','linewidth',2);
    dnG = size(nameFig,2);

    i = reco{1};
    visualisationShared(PSI_z*mu_w{i}, PSI_z*1.96*sqrt(diag(sigma_w{i} )), sum(nbDof), z,  i, 'b', nameFig);
    nameFig = visualisation(PSI_z*mu_w{i}, sum(nbDof), z, i, 'b', nameFig);
    prevG = size(nameFig,2);
    visualisationShared(PSI_z*mu_new, PSI_z*1.96*sqrt(diag(sigma_new)), sum(nbDof), z,  i, 'g', nameFig);
    nameFig = visualisation(PSI_z*mu_new, sum(nbDof), z, i,'g', nameFig);
    newG = size(nameFig,2);
    switch(vff)
        case 1
            title(['Inference of trajectory ', num2str(trial)], 'fontsize',30 )
            ylabel('X (m)', 'fontsize', 24);
        case 2
            ylabel('Y (m)', 'fontsize', 24);
        case 3
            ylabel('Z (m)', 'fontsize', 24);
        case 4
            ylabel ('Fx (N)', 'fontsize', 24);
        case 5 
             ylabel ('Fy (N)', 'fontsize', 24);
             xlabel('Iterations', 'fontsize', 24);
        case 6 
             ylabel ('Fz (N)', 'fontsize', 24);
             xlabel('Iterations', 'fontsize', 24);
    end
end
legend(nameFig(1,[dtG, dnG, prevG, newG]),'desired trajectory', 'data known','learnt distribution', 'new distribution' );
