%here we compute the error between courb and so on
cpt=1;
i=1
figAnalysis = figure;


   for type=1:nbDofTot
  %   subplot(nb +1,(nbDofTot-nb)+1,type) 
        subplot(3,1,1);
        cpt=1;
        a = size(figAnalysis,2);
        
        for val=1:size(newmu,2)
            %tallActu = size(nf3D,2);
            k=1+  floor(totalTimeTrial(trial) /size(newmu,2))*(val-1);
     %       figAnalysis((a + (val-1)*2 + 1))=  plot(val*(z /size(newmu,2)), y_trial_Tot{i}(val*(totalTimeTrial(trial) /size(newmu,2)) + (type-1)*totalTimeTrial(i)), '*b');%, 'Color' , [0, val/11, 0] );
            figAnalysis = visualisation(PSI_z*newmu{val}, nbDofTot, z,type, [0, val/size(newmu,2), 0], figAnalysis);
        end
        figAnalysis = visualisationShared(PSI_z*newmu{size(newmu,2)}, PSI_z*1.96*sqrt(diag(newSigma{size(newmu,2)})), nbDofTot, z,type, 'g', figAnalysis);
        b = size(figAnalysis,2);
       figAnalysis = visualisation2(y_trial_Tot{i}, nbDofTot, totalTimeTrial(i),type, '--magenta', realAlpha, figAnalysis);hold on;
     %  figAnalysis = visualisation2(y_trial_Tot{i}, nbDofTot, totalTimeTrial(i),type, '-.r', mu_alpha(i), figAnalysis);hold on;
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
 legend(figAnalysis([(a+2) b (b+1) ]), 'first update', 'last update', 'trial');%,% 'trial with the supposed velocity'); (b+2)
    annotation('textbox', [.3 .9 .1 .1], 'String', ['acc=', num2str(accuracy), '. nbFct=', num2str(nbFunctions(1))]); 




    for i=1:nbDofTot
        for j=1:totalTimeTrial
            val(i,j) = y_trial_Tot{1}(totalTimeTrial*(i-1)+j);
        end
    end
   % plot(val(1,:),'.-r');hold on;
    subplot(3,1,2);
clear logLike;
tee = size(figAnalysis,2);
for test=1: floor(totalTimeTrial /size(newmu,2)): totalTimeTrial(1)
    testC = PSI_z*newmu{cpt};
    testS = PSI_z*newSigma{cpt}*PSI_z';
    %loglike
    d =  det(testS+ 0.001*eye(size(testS,2)));
    a(cpt) =(val(1,:) - testC')*inv(testS)*(val(1,:) - testC')';
    vall(cpt)= sum((val(1,:) - testC'))*sum((val(1,:) - testC'))';
    %a=  trace(inv(testS)*(val(1,:) - testC')'*(val(1,:) - testC'))
    %a*inv(a) 
    logLike(cpt) = -(100/2)*log(2*pi)  -0.5*a(cpt) -0.5*log(d)
    %logLikelihood(y_trial_Tot{1},testC,testS) 
    figAnalysis(size(figAnalysis,2)+1) = plot(abs(testC'-val(1,:)), 'Color', [0, cpt/size(newmu,2), 0]);hold on;
    plot(100*(cpt /size(newmu,2) ), 0.05, '+r');
    cpt= cpt+1;
end

xlabel('Iteration');
ylabel('distance (m): y_d - y_{infered}');
title('Error between the real trajectory and the infered ones');
legend(figAnalysis([(tee+1), size(figAnalysis,2)]), 'first inference', 'last inference')


subplot(3,1,3);
plot(logLike, '-o');hold on;
title('log likelihood according to time');
xlabel('iterations')
ylabel('log likelihood')
%xaxis([0 100]);