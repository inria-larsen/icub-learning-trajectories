f = figure;
    set(0,'DefaultAxesFontSize',18)
f(1:5) = plot(PSI_z, 'k');hold on;
basis = size(f,2);
f(size(f,2) +1) = plot(y{1}{1}(1:(totalTime(1,1)/ z):totalTime(1,1),1), '+');hold on;
data = size(f,2);
f(size(f,2) +1) = plot(PSI_z*w{1}(1,:)','r');hold on;
infered = size(f,2);

for i=1:5
   f(size(f,2) +1) = plot(PSI_z(:,i)*w{1}(1,i), '--g');hold on
end

legend(f([basis data infered size(f,2)]), 'basis functions', 'real data', 'computed trajectory', 'w(i)phi(i)')
 xlabel('iterations')
  ylabel('x')
   title('Example of the computation for one trajectory')

% 
% f = figure;
% set(0,'DefaultAxesFontSize',18)
% for i=1:30
%     f(size(f,2) +1) = plot(y{1}{i}(1:(totalTime(1,i)/ z):totalTime(1,i),1), 'g');hold on;
% end
% data = size(f,2)
% f =  visualisationShared(PSI_z*mu_w{size(mu_w,2)}, PSI_z*1.96*sqrt(diag(sigma_w{size(mu_w,2)})), 1, z,1, 'r', f);
% f(size(f,2) +1) = plot(PSI_z*mu_w{size(mu_w,2)},'r');
% infered = size(f,2)
% 
% for i=1:5
%     f(size(f,2) +1) = plot(PSI_z(:,i)*w{1}(1,i), '--k');hold on
% end
% 
% legend(f([data infered size(f,2)]), 'real data', 'computed distribution', 'w(i)phi(i)')
% xlabel('iterations')
% ylabel('x')
% title('Example of a computed distribution of a trajectory')