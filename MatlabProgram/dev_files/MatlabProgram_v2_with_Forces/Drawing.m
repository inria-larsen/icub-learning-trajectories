% Plot the trials we will use to verify the learning

close all
clear nam*

%%%PLOT DISTRIBUTION
    nf3D = figure(44);
%    % Plot2D the learned distribution 
%    type=1;
%    nf3D = visualisationShared(PSI_z*mu_w{i},PSI_z*1.96*sqrt(diag(sigma_w{i})), nbDofTot, z, type, 'r', nf3D); hold on;
%    nf3D = visualisation(PSI_z*mu_w{i}, nbDofTot, z, i,'r', nf3D);
%    a=size(nf3D,2)
%    for j=1:var(i)
%     	nf3D = visualisation2(y{i}{j},nbDofTot, totalTime(i,j), type, 'b',alpha2{i}(j), nf3D);
%     end 
%     title(['learned distribution for the' num2str(type) 'th DOF of the' num2str(i) 'th trajectory' ]);%strcat('learned distribution for the', atostr(i), ' trajectory');
%     xlabel(['x'])
%     legend(nf3D([a (size(nf3D,2)-1)]),'Distribution Learned','Data');
%        



%%%%%%%%%%%%%%END PLOT UPDATE

   i=1
     %nb = mod(nbDofTot,3);
   % Plot2D the learned distribution 
   for type=1:nbDofTot
     subplot(3,2,type);
      %subplot(nb,(nbDofTot-nb)+1,type) 
    
       nf3D = visualisationShared(PSI_z*mu_w{i},PSI_z*1.96*sqrt(diag(sigma_w{i})), nbDofTot, z, type, 'r', nf3D); hold on;
                       title(['Learned distribution with acc=', num2str(accuracy), ', nbFunctions =', num2str(nbFunctions(1))]);

       nf3D = visualisation(PSI_z*mu_w{i}, nbDofTot, z, type,'r', nf3D);
        a=size(nf3D,2);
       
%   
%         nf3D = visualisationShared(PSI_z*newmu{11}, PSI_z*1.96*sqrt(diag(newSigma{11})), nbDofTot, z,type, 'g', nf3D);
%         nf3D = visualisation(PSI_z*newmu{11}, 6, z,type, [0, 1, 0], nf3D);
%      
        b = size(nf3D,2);
       
       for j=1:var(i)
            nf3D = visualisation2(y{i}{j},nbDofTot, totalTime(i,j), type, ':b',alpha2{i}(j), nf3D);
       end 
       c = size(nf3D,2);
      % nf3D = visualisation2(y_trial_Tot{i}, nbDofTot, totalTimeTrial(i),type, '--magenta', realAlpha, nf3D);hold on;
       if (type ==1) 
                  nm= 'x';
       elseif(type ==2) nm = 'y';
       elseif(type == 3) nm= 'z';
       elseif(type == 4) nm = 'fx';
       elseif(type == 5) nm = 'fy';
       else nm = 'fz';
       end
        xlabel('iterations');
        ylabel(nm);
   end
    legend(nf3D([ a b c]),'Distribution Learned', 'update', 'data');
     annotation('textbox', [.3 .9 .1 .1], 'String', ['NbFct=', num2str(nbFunctions(1))]);
    
     
     %%update acc to time
    clear nf3D;
   nf3D = figure(3);
   i=1;
   % Plot2D the learned distribution 
   for type=1:nbDofTot
  %   subplot(nb +1,(nbDofTot-nb)+1,type) 
        subplot(3,2,type);
        cpt=1;
        a = size(nf3D,2);
        val = 1
        %for val=1:2:size(newmu,2)
            %tallActu = size(nf3D,2);
            k=1+  floor(min(z,totalTimeTrial(trial)) /size(newmu,2))*(val-1);
            nf3D((a + (val-1)*2 + 1))=  plot(val*(z /size(newmu,2)), y_trial_Tot{i}(val*floor(totalTimeTrial(trial) /size(newmu,2)) + (type-1)*totalTimeTrial(i)), '*b');hold on;%, 'Color' , [0, val/11, 0] );
            nf3D = visualisation(PSI_z*newmu{val}, nbDofTot, z,type, [0, val/size(newmu,2), 0], nf3D);
       % end
       val = size(newmu,2);
            k=1+  floor(min(z,totalTimeTrial(trial)) /size(newmu,2))*(val-1);
            nf3D((a + (val-1)*2 + 1))=  plot(val*(z /size(newmu,2)), y_trial_Tot{i}(val*floor(totalTimeTrial(trial) /size(newmu,2)) + (type-1)*totalTimeTrial(i)), '*b');hold on;%, 'Color' , [0, val/11, 0] );
            nf3D = visualisation(PSI_z*newmu{val}, nbDofTot, z,type, [0, val/size(newmu,2), 0], nf3D);
       
        %nf3D = visualisationShared(PSI_z*newmu{size(newmu,2)}, PSI_z*1.96*sqrt(diag(newSigma{size(newmu,2)})), nbDofTot, z,type, 'g', nf3D);
        b = size(nf3D,2);
       nf3D = visualisation2(y_trial_Tot{i}, nbDofTot, totalTimeTrial(i),type, '--magenta', realAlpha, nf3D);hold on;
       nf3D = visualisation2(y_trial_Tot{i}, nbDofTot, meanTime(i) ,type, '-.r', mu_alpha(i), nf3D);hold on;
        if (type ==1) nm= 'x';
           elseif(type ==2) nm = 'y';
           elseif(type == 3) nm= 'z';
           elseif(type == 4) nm = 'fx';
           elseif(type == 5) nm = 'fy';
           else nm = 'fz';
       end
       ylabel(nm);
       xlabel('iterations');
   end
 legend(nf3D([(a+2) b (b+1) (b+2)]), 'first update', 'last update', 'trial', 'trial with the supposed velocity'); 
    annotation('textbox', [.3 .9 .1 .1], 'String', ['acc=', num2str(accuracy), '. nbFct=', num2str(nbFunctions(1))]);    
    % subplot(nb +1,(nbDofTot-nb)+1,nbDofTot+1) 
    
    figure;
    %subplot(3,1,3);
   
    for t=1:size(newmu,2)
      scatter(t*(z /size(newmu,2)), prob(i,t), '+');hold on;
    end
    Xax =[(z /size(newmu,2)):(z /size(newmu,2)):100];
    plot(Xax, prob(i,:));hold on;
    title('log likelihood according to time');
    xlabel('iterations')
    ylabel('log likelihood')
% % plot 3D trials
% % nam4 = figure;
% % for k =1 : nbKindOfTraj
% %     for i=1:var(k)
% %       nam4 = visualisation3D(y_trial_Tot{k}, nbDofTot, totalTimeTrial(k), 2, nbDof, [(nbKindOfTraj-k)/nbKindOfTraj, k/nbKindOfTraj,0.5*(k/nbKindOfTraj)], nam4);hold on;
% %     end
% % end
% % title('Trial trajectories');
% % xlabel('Fx')
% % ylabel('Fy')
% % zlabel('Fz')
% 
% 
%   %close all
% %   %clear namf
% %   
% %   nbData = reco{12}(3);
% %   namf = figure;    
% %     namf = visualisationShared(PSI_z*newmu{11}, PSI_z*1.96*sqrt((diag(newSigma{11}))), 6, z,1, 'g', namf);
% %       namf = visualisation(PSI_z*newmu{11}, 6, z,1, 'g', namf);
% %   b = size(namf,2);
% % 
% % 
% %   for t=1:nbData
% %            namf(1+ size(namf,2)) = scatter(t*realAlpha, y_trial_nbData((i-1)*nbData + t), '.b'); hold on;       
% %     end
% %     namf = visualisation2(y_trial_Tot{1}, 6, totalTimeTrial(1), i, ':b', realAlpha, namf);
% % 
% %     
% %     legend(namf([a  b (b+1) (b + nbData + 1) ]),'prev Learned distribution', 'new distribution', 'Data known', 'Data we should have');
% %     title(['Position recognition of the trajectory ']);
% %     xlabel('Iterations');
% %     ylabel('x');
% 
% 
% 
% %      % Plot3D the learned distribution 
% %    for type=1:size(nbDof,2)
% %     
% %          clear nf3D
% %        nf3D = figure;    
% %        for j=1:3:var(i)
% %            nf3D = visualisation3D(y{i}{j},nbDofTot, totalTime(i,j), type, nbDof, ':b', nf3D);
% %        end
% %        nf3D = visualisation3D(PSI_z*mu_w{i}, nbDofTot, z, type, nbDof, '.r', nf3D); hold on;
% %         nf3D = visualisation3D(PSI_z*mu_w{i}, nbDofTot, z, type, nbDof, 'r', nf3D); hold on;
% %         nf3D = visualisation3D(PSI_z*mu_w{i}, nbDofTot, z, type, nbDof, 'or', nf3D); hold on;
% % 
% %        nf3D = visualisation3D(PSI_z*(mu_w{i} + 1.96*sqrt(diag(sigma_w{i}))), nbDofTot, z,type, nbDof, '-.g', nf3D);
% %        nf3D = visualisation3D(PSI_z*(mu_w{i}- 1.96*sqrt(diag(sigma_w{i}))),nbDofTot, z,type, nbDof, '-.g', nf3D);
% %         title(['learned forces distribution' ]);%strcat('learned distribution for the', atostr(i), ' trajectory');
% %         
% %         xlabel('Fx')
% %         ylabel('Fy')
% %         zlabel('Fz')
% %         legend(nf3D([2 (size(nf3D,2)-3)]), 'data', 'distribution')   
% %    end
% 
% 
% 
% % nbData = reco{cpt}(1);
% % nameFig = figure;
% %   %draw the infered movement
% %       for t=1:nbData
% %            nameFig(t) = scatter(t*realAlpha, y_trial_part((3-1)*nbData + t), '.b'); hold on;       
% %     end
% %     nameFig = visualisation2(y_trial_Tot, 6, totalTimeTrial, i, ':b', realAlpha, nameFig);
% % 
% %     nameFig = visualisationShared(PSI_z*newmu{1}, PSI_z*1.96*sqrt((diag(newSigma{1})), 6, z,1, 'g', nameFig);
% %     nameFig = visualisation(PSI_z*mu_w{i},PSI_z*1.96*sqrt(diag(sigma_w{i})) , 6, z, i,'r', nameFig);
% %     legend(nameFig([1 (nbData+1) (nbData +2) (nbData +3) ]),'Data known', 'Data we should have', 'Data deducted', 'prev Learned distribution');
% %     title(['Position recognition of the trajectory ']);
% %     xlabel('Iterations');
% %     ylabel('x');
% 
% % figure;
% % for t=1:size(newmu,2)
% %     scatter(t, prob(i,t), '+');hold on;
% % end
% % title('log likelihood according to time');
% % xlabel('iterations')
% % ylabel('log likelihood')
% 
i=1;
fig = figure;
filename = 'testGif';
           for type=1:nbDofTot
               subplot(3,2,type);
           %     fig = visualisation2(y_trial_Tot{i}, nbDofTot, totalTimeTrial(i),type, '--magenta', realAlpha, fig);hold on;
                fig = visualisation2(y_trial_Tot{i}, nbDofTot, meanTime(i),type, '-.r', mu_alpha(i), fig);hold on;
           end
                               drawnow
            frame = getframe(1);

             im = frame2im(frame);
             [imind,cm] = rgb2ind(im,256);
                 imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
        for val=1:size(newmu,2)
            for type=1:nbDofTot
                subplot(3,2,type);
                cpt=1;
            %tallActu = size(nf3D,2);
            k=1+  floor(totalTimeTrial(trial) /size(newmu,2))*(val-1);
            plot(val*(z /size(newmu,2)), y_trial_Tot{1}(val*floor(totalTimeTrial(trial) /size(newmu,2)) + (type-1)*totalTimeTrial(trial)), '*', 'Color' , [0, val/11, 0] );
            visualisation(PSI_z*newmu{val}, nbDofTot, z,type, [0, val/size(newmu,2), 0], fig);
           if (type ==1) nm= 'x';
           elseif(type ==2) nm = 'y';
           elseif(type == 3) nm= 'z';
           elseif(type == 4) nm = 'fx';
           elseif(type == 5) nm = 'fy';
           else nm = 'fz';
       end
      ylabel(nm);
      xlabel('iterations');
        end
                    drawnow
            frame = getframe(1);

             im = frame2im(frame);
             [imind,cm] = rgb2ind(im,256);
             imwrite(imind,cm,filename,'gif','WriteMode','append');
        fig = visualisationShared(PSI_z*newmu{size(newmu,2)}, PSI_z*1.96*sqrt(diag(newSigma{size(newmu,2)})), nbDofTot, z,type, 'g', fig);
        b = size(fig,2);
    annotation('textbox', [.3 .9 .1 .1], 'String', ['acc=', num2str(accuracy), '. nbFct=', num2str(nbFunctions(1))]);   
 
        drawnow
            frame = getframe(1);

             im = frame2im(frame);
             [imind,cm] = rgb2ind(im,256);
             
                 imwrite(imind,cm,filename,'gif','WriteMode','append');

   end

  
