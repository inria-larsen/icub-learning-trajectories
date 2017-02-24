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
    min_alpha_i(i) = min(alpha{i});
    max_alpha_i(i) = max(alpha{i});
end
min_alpha = min(min_alpha_i(i));
max_alpha = max(max_alpha_i(i));

PSI_z = computeBasisFunction (z,nbFunctions,1, z,h);

%w computation for each trials
%figure;
for i = 1 : nbKindOfTraj
    test = zeros(var(i),nbFunctions(1)*nbDof(1)+nbFunctions(2)*nbDof(2));
    
    
    for j = 1 : var(i)
        %resolve a little bug
        sizeY  = size(y{i}{j},1);
        if(sizeY ~= size(PSI{i}{j},1))
            y{i}{j} = y{i}{j}(1:sizeY-(nbDof(1)+nbDof(2)));
            totalTime(i,j) = totalTime(i,j) -6;
            alpha{i}(j) = z /totalTime(i,j);
        end

        w{i}(j,:) = (PSI{i}{j}'*PSI{i}{j}+1e-6*eye(nbFunctions(1)*nbDof(1)+nbFunctions(2)*nbDof(2))) \ PSI{i}{j}' * y{i}{j};        
        test(j,:) =w{i}(j,:); 
    end
        
    %computation of the w distribution     
    mu_w{i} = mean(test)';
    sigma_w{i} = cov(test);

   % Plot the learned distribution for the x cartesian position of the
   % trajectory i
   nf = figure;         
   nf = visualisation(PSI_z*mu_w{i}, 6, z, 1, 'r', nf);
   nf = visualisation(PSI_z*(mu_w{i} + 1.96*sqrt(diag(sigma_w{i}))), 6, z, 1,'-.r', nf);
   nf = visualisation(PSI_z*(mu_w{i}- 1.96*sqrt(diag(sigma_w{i}))), 6, z, 1, '-.r', nf);
   for j=1:var(i)
       nf = visualisation2(y{i}{j}, 6, totalTime(i,j),1, 'b', alpha{i}(j) , nf);hold on;
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
   
   
     % Plot3D the learned distribution for the forces position of the
   % trajectory i
   nf3D = figure;         
   nf3D = visualisation3D(PSI_z*mu_w{i}, 6, z, 1, 'r', nf3D); hold on;
   nf3D = visualisation3D(PSI_z*(mu_w{i} + 1.96*sqrt(diag(sigma_w{i}))), 6, z,1, '-.g', nf3D);
   nf3D = visualisation3D(PSI_z*(mu_w{i}- 1.96*sqrt(diag(sigma_w{i}))), 6, z,1, '-.g', nf3D);
   for j=1:var(i)
       nf3D = visualisation3D2(y{i}{j}, 6, totalTime(i,j), 1, 'b', alpha{i}(j) , nf3D);
   end 
    
    if(i==1)
        title('learned distribution for the forces of the trajectory right');%strcat('learned distribution for the', atostr(i), ' trajectory');
    elseif(i==2)
        
        title('learned distribution for the forces of the trajectory ahead');   
    else
        title('learned distribution for the forces of the trajectory top');
    end
    xlabel('Fx (N)')
    ylabel('Fy (N)')
    zlabel('Fz (N)')
   legend(nf3D([2 4 5]),'Distribution Learned','Standard deviation','Data'); 
   
     % Plot3D the learned distribution for the cartesian position of the
   % trajectory i
   nf3D = figure;         
   nf3D = visualisation3D(PSI_z*mu_w{i}, 6, z, 0, 'r', nf3D); hold on;
   nf3D = visualisation3D(PSI_z*(mu_w{i} + 1.96*sqrt(diag(sigma_w{i}))), 6, z,0, '-.g', nf3D);
   nf3D = visualisation3D(PSI_z*(mu_w{i}- 1.96*sqrt(diag(sigma_w{i}))), 6, z,0, '-.g', nf3D);
   for j=1:var(i)
       nf3D = visualisation3D2(y{i}{j}, 6, totalTime(i,j), 0, 'b', alpha{i}(j) , nf3D);
   end 
    
    if(i==1)
        title('learned distribution for the position of the trajectory right');%strcat('learned distribution for the', atostr(i), ' trajectory');
    elseif(i==2)
        
        title('learned distribution for the position of the trajectory ahead');   
    else
        title('learned distribution for the position of the trajectory top');
    end
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
   legend(nf3D([2 4 5]),'Distribution Learned','Standard deviation','Data'); 
   
end
    clear i j test w PSI min_alpha_i max_alpha_i;
