 list = {'x','y','z','fx','fy','fz'};
   	

for i=1:nbKindOfTraj
    fig = figure(i+nbDof(1));
 
    for vff=1:nbDofTot
            subplot(nbDof(1),size(nbDof,2),vff);
%             a = PSI_z*mu_w{i};
%             b = PSI_z*1.96*sqrt(diag(sigma_w{i}));% a +
%             plot(a(1+100*(vff-1):100*(vff)),'g');hold on;
%              plot(a(1+100*(vff-1):100*(vff)) - b(1+100*(vff-1):100*(vff)), '.-g');hold on;
%              plot(a(1+100*(vff-1):100*(vff)) + b(1+100*(vff-1):100*(vff)), '.-g');hold on;
            %fig(size(fig,2) + 1) =
            fig = visualisationShared(PSI_z*mu_w{i}, PSI_z*1.96*sqrt(diag(sigma_w{i} )), nbDofTot, z,  vff, 'b', fig);
            fig = visualisation(PSI_z*mu_w{i}, nbDofTot, z,  vff, 'b', fig);
            disG = size(fig,2);
            for j = 1 : var(i)
                %a = w{i}(j,:)*PSI_z';
                %plot(a((vff-1)*100 + 1:(vff)*100) , ':k');hold on;
              fig(size(fig,2) + 1) =  plot(y{i}{j}(1 + totalTime(i,j)*(vff-1) : totalTime(i,j)/100 : totalTime(i,j)*vff), ':k','linewidth',2);hold on;
            end
            datG = size(fig,2);
            %un des tests
            %plot(y_trial_Tot{i}(totalTimeTrial(i)*(vff-1) + 1 : (totalTimeTrial(i)/100): totalTimeTrial(i)*vff),'+g');hold on;
        %visualisationShared(a((vff-1)*100 + 1:(vff)*100) , b((vff-1)*100 + 1:(vff)*100 ), 100, 1, z, 1, fig);hold on;
       % plot(b((vff-1)*100 + 1:(vff)*100) ,'--', 'Color', [i/nbKindOfTraj, 0, 0]);hold on;
       % b = a - PSI_z*1.96*sqrt(diag(sigma_w{i}));
       % plot(b((vff-1)*100 + 1:(vff)*100) ,'--', 'Color', [i/nbKindOfTraj, 0, 0]);hold on;
         xlabel('iterations', 'fontsize', 24);
         ylabel(list{vff}, 'fontsize', 24);
         if(vff==1)    
             title(['trajecories type ', num2str(i)], 'fontsize', 30);
         end
    end
        
    
        legend(fig([disG,datG]), 'distribution learnt (mean & standart deviation)','observed data');
end