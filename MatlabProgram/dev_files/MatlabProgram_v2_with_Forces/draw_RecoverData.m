list = {'x','y','z','fx','fy','fz'};

for k =1:nbKindOfTraj
    fig22 = figure;
    for l=1:nbDofTot  
        subplot(nbDof(1),size(nbDof,2),l);
        for i=1:var(k)      
            fig22 = visualisation(y{k}{i},nbDofTot,totalTime(k,i), l, ':k',fig22);hold on;
        end
        fig22=visualisation(y_trial{k}, nbDofTot,nbData,l, '.b',fig22);hold on;
         fig22=visualisation(y_trial_Tot{k}, nbDofTot, totalTimeTrial(k),l, 'b',fig22);hold on;
         xlabel('iterations', 'fontsize', 24);
         ylabel(list{l}, 'fontsize', 24);
         if(l==1)    
             title(['trajecories type ', num2str(k)], 'fontsize', 30);
         end
    end
end


