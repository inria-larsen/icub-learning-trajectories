% In this part, we compute the distribution for each kind of trajectories.

for i=1:nbKindOfTraj
   
    for j = 1:var(i)
        %we compute the phasis
        alpha{i}(j) = z / totalTime(i,j);
        %we compute the corresponding PSI matrix
        PSI{i}{j} = computeBasisFunction (z,nbFunctions,alpha{i}(j), totalTime(i,j),h);
    end
    mu_alpha(i) = mean(alpha{i});
    sigma_alpha(i) = cov(alpha{i});
end

PSI_z = computeBasisFunction (z,nbFunctions,1, z,h);

%w computation for each trials
%figure;
for i = 1 : nbKindOfTraj
    test = zeros(var(i),nbFunctions(1)*nbDof(1)+nbFunctions(2)*nbDof(2));
    for j = 1 : var(i)
		w{i}(j,:) = (PSI{i}{j}'*PSI{i}{j}+1e-6*eye(nbFunctions(1)*nbDof(1)+nbFunctions(2)*nbDof(2))) \ PSI{i}{j}' * y{i}{j};
		test(j,:) =w{i}(j,:); 
    end
        
    %computation of the w distribution     
    mu_w{i} = mean(test)';
    sigma_w{i} = cov(test);

   % Plot the learned distribution for the x cartesian position of the
   % trajectory i
   nf = figure;         
   nf = visualisation(PSI_z*mu_w{i}, sum(nbDof), z, i, 'r', nf);
   nf = visualisation(PSI_z*(mu_w{i} + 1.96*sqrt(diag(sigma_w{i}))),sum(nbDof), z, i, '-.r', nf);
   nf = visualisation(PSI_z*(mu_w{i}- 1.96*sqrt(diag(sigma_w{i}))), sum(nbDof), z, i,'-.r', nf);
   for j=1:var(i)
       nf = visualisation2(y{i}{j}, sum(nbDof), totalTime(i,j), i, 'b', alpha{i}(j) , nf);hold on;
   end 
    
    if(i==1)
        title('learned distribution for the x position of the trajectory right');%strcat('learned distribution for the', atostr(i), ' trajectory');
    elseif(i==2)
        
        title('learned distribution for the x position of the trajectory ahead');   
    else
        title('learned distribution for the x position of the trajectory top');
    end
    xlabel('Iteration')
    ylabel('x (m)')
    legend(nf([2 4 5]),'Distribution Learned','Standard deviation','Data'); 
end
    clear i j test w PSI;
