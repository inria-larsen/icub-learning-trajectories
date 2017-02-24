%%%%%%%%%%%%%%%%%%%%%%%%%%%% Recover the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nbKindOfTraj = 3;
name{1} = 'Data/dataAhead.txt';
name{2} = 'Data/dataTop.txt';
name{3} = 'Data/dataRight.txt';
for i=1:nbKindOfTraj
% %we open the files
f(i) = fopen(name{i}, 'r');

% we scan the files line per line
X{i} = textscan(f(i), '%s', 'delimiter', '\n');

%we close the files
fclose(f(i));

%we put data in vectors as numbers
X{i} = cellfun(@str2num, X{i}{1}, 'UniformOutput', false);

%we compute the number of block (delimited by a carriage return) and we put
%data in the vector B.
numBlock = 1;
numLine = 0;
for n = 1:size(X{i}, 1)
 
    if isempty(X{i}{n})
        numBlock = numBlock + 1;
        numLine = 0;
    else
        numLine = numLine+1;
        data{i}{numBlock}(numLine,:) = X{i}{n}; 
    end
 
end

end
clear f numBlock numLine n X ans
% load('dataTotal2.mat');
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create the y vector  and other usefull variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

a = 1;
v = [1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10];
namf = figure;
affich=1; %what we want to plot  (x y z fx fy fz))

% % We keep the last sample to try to recognize the movement
for k=1:nbKindOfTraj
    var(k) = size(data{k},2) - 1;
    % totalTime will be the number of time of each trajectories, to begin with it is the same as the totalTime
    % y is the input vector of data
    for i=1:var(k)
        y{k}{i} = [ data{k}{i}(:,1) ; data{k}{i}(:,2); data{k}{i}(:,3); filter(v,a,data{k}{i}(:,4)) ; filter(v,a,data{k}{i}(:,5)); filter(v,a,data{k}{i}(:,6))];
        totalTime(k,i) = size(data{k}{i},1);
        namf = visualisation(y{k}{i}, 6, totalTime(k,i), affich, [(nbKindOfTraj-k)/nbKindOfTraj, k/nbKindOfTraj, 0], namf); %we draw only the first DOF because too much data  
    end 
end
% for i=1:var(2)
%     totalTime(2,i) = size(Ba{i},1);
%     y{2}{i} = [ Ba{i}(:,1) ; Ba{i}(:,2); Ba{i}(:,3) ;  filter(v,a,Ba{i}(:,4)) ; filter(v,a,Ba{i}(:,5)); filter(v,a,Ba{i}(:,6))];
%   %  namf = visualisation(y{2}{i}, 6, totalTime(2,i),affich,  [0,0,1], namf );
%     
% end
% 
% for i=1:var(3)
%     totalTime(3,i) = size(Bt{i},1);
%     y{3}{i} = [ Bt{i}(:,1) ; Bt{i}(:,2); Bt{i}(:,3) ;  filter(v,a,Bt{i}(:,4)) ;  filter(v,a,Bt{i}(:,5));  filter(v,a,Bt{i}(:,6))];
%   %  namf = visualisation(y{3}{i}, 6, totalTime(3,i),affich,  [0,1,0], namf);
%     
% end
% 
% %Information about the previous plot that represent each sample of each
% %trajectory
% % if(affich ==1) namePlot = ('x');
% % elseif(affich ==2) namePlot = ('y');
% % elseif(affich ==3) namePlot = ('z');
% % elseif(affich ==4) namePlot = ('fx');
% % elseif(affich ==5) namePlot = ('fy');
% % else namePlot = ('fz');
% % end
% % title('Trajectories used to learn');
% % xlabel('time');
% % ylabel(namePlot);
% % legend(namf([2 (var(1)+2) (var(1) + var(2) + 2)]),'Right','Ahead','Top');
% 
% 
% 
% %TRIAL DATA
% %we keep the last sample of each trajectory to see if we recognize
% %correctly the last movement
% % i = size(Br,2);
% % totalTimeTrial(1) = size(Br{i},1); %if we continue the movement to the end
% % y_trial_Tot{1} =  [ Br{i}(:,1) ; Br{i}(:,2); Br{i}(:,3); filter(v,a,Br{i}(:,4)) ; filter(v,a,Br{i}(:,5)); filter(v,a,Br{i}(:,6))];
% % %the first data of the trajectory from where we would have to recognize the movement
% % y_trial{1} = [ Br{i}(1:nbData,1) ; Br{i}(1:nbData,2); Br{i}(1:nbData,3) ; filter(v,a,Br{i}(1:nbData,4)) ;filter(v,a, Br{i}(1:nbData,5)); filter(v,a,Br{i}(1:nbData,6))]; 
% % 
% % i = size(Ba,2);
% % totalTimeTrial(2) = size(Ba{i},1);
% % y_trial_Tot{2} =  [ Ba{i}(:,1) ; Ba{i}(:,2); Ba{i}(:,3) ;  filter(v,a,Ba{i}(:,4)) ; filter(v,a,Ba{i}(:,5)); filter(v,a,Ba{i}(:,6))];
% % y_trial{2} = [ Ba{i}(1:nbData,1) ; Ba{i}(1:nbData,2); Ba{i}(1:nbData,3) ; filter(v,a,Ba{i}(1:nbData,4)) ; filter(v,a,Ba{i}(1:nbData,5)); filter(v,a,Ba{i}(1:nbData,6))];
% % 
% % i = size(Bt,2);
% % totalTimeTrial(3) = size(Bt{i},1);
% % y_trial_Tot{3} = [ Bt{i}(:,1) ; Bt{i}(:,2); Bt{i}(:,3) ;  filter(v,a,Bt{i}(:,4)) ;  filter(v,a,Bt{i}(:,5));  filter(v,a,Bt{i}(:,6))];
% % y_trial{3} = [ Bt{i}(1:nbData,1) ; Bt{i}(1:nbData,2); Bt{i}(1:nbData,3) ; filter(v,a,Bt{i}(1:nbData,4)) ; filter(v,a,Bt{i}(1:nbData,5)) ; filter(v,a,Bt{i}(1:nbData,6))];
% % 
% % 
% % for j=1:nbKindOfTraj 
% %     for t=1:nbData
% %          for i=1:nbDof(1) %+ nbDof(2)
% %              yaff{j}((nbData*(i-1) + t),1) = y_trial_Tot{j}((totalTimeTrial(j)*(i-1) + t),1); %here we have only the coordinate data
% %          end
% %     end
% % end
% 
% % % Plot the trials we will use to verify the learning
% % name = figure;
% % name = visualisation(y_trial_Tot{1}, 6, totalTimeTrial(1),affich, [1, 0, 0], name);hold on;
% % name = visualisation(y_trial_Tot{2},6, totalTimeTrial(2),affich, [0, 0, 1], name);
% % name = visualisation(y_trial_Tot{3},6, totalTimeTrial(3),affich, [0, 1, 0], name);
% % name = visualisation(y_trial{1}, 6, nbData,affich, '.r', name);
% % name = visualisation(y_trial{2}, 6, nbData,affich, '.b', name);
% % name = visualisation(y_trial{3}, 6, nbData,affich, '.g', name);
% % title('Trial we keep to try the recognition step (one for each kind of trajectory)');
% % xlabel('Iterations');
% % ylabel(namePlot);
% % legend('Right','Ahead','Top');
% % 
% % 
% % %%%%%%%%%%%%%%%%%%%3D REPRESENTATION
% 
% % %plot 3D data
% % nam2 = figure;
% % for i=1:var(1)
% %     nam2 = visualisation3D(y{1}{i}, 6, totalTime(1,i), 0, [1, 0, 0], nam2);
% % end
% % for i=1:var(2)
% %     nam2 = visualisation3D(y{2}{i}, 6, totalTime(2,i), 0,[0, 0, 1], nam2);
% % end
% % for i=1:var(3)
% %     nam2 = visualisation3D(y{3}{i}, 6, totalTime(3,i),0,[0, 1, 0], nam2);
% % end
% % 
% % nam2 = visualisation3D(y_trial_Tot{1}, 6, totalTimeTrial(1), 0,  [0.5, 0.5, 0.5], nam2);
% % nam2 = visualisation3D(y_trial_Tot{2}, 6, totalTimeTrial(2), 0, [0.5, 0.5, 0.5], nam2);
% % nam2 = visualisation3D(y_trial_Tot{3}, 6, totalTimeTrial(3),0, [0.5, 0.5, 0.5], nam2);
% % nam2 = visualisation3D(y_trial{1}, 3, nbData, 0, '.b', nam2);
% % nam2 = visualisation3D(y_trial{2}, 3, nbData, 0, '.b', nam2);
% % nam2 = visualisation3D(y_trial{3}, 3, nbData, 0,'.b', nam2);
% % 
% % title('Trajectories used to learn')
% % xlabel('x')
% % ylabel('y')
% % zlabel('z')
% % legend(nam2([2 (var(1)+2) (var(1) + var(2) + 2) (var(1) + var(2) + var(3) + 2)]),'Right','Ahead','Top', 'Trials');
% 
% %plot 3D data forces
% nam3 = figure;
% for i=1:var(1)
%     nam3 = visualisation3D(y{1}{i}, 6, totalTime(1,i), 1, [1, 0, 0], nam3);
% end
% for i=1:var(2)
%     nam3 = visualisation3D(y{2}{i}, 6, totalTime(2,i), 1,[0, 0, 1], nam3);
% end
% for i=1:var(3)
%     nam3 = visualisation3D(y{3}{i}, 6, totalTime(3,i),1,[0, 1, 0], nam3);
% end
% title('Forces used to learn')
% xlabel('Fx')
% ylabel('Fy')
% zlabel('Fz')
% legend(nam3([2 (var(1)+2) (var(1) + var(2) + 2)]),'Right','Ahead','Top');
% 
% % %plot 3D trials
% % nam4 = figure;
% % nam4 = visualisation3D(y_trial_Tot{1}, 6, totalTimeTrial(1), 0, [1, 0, 0], nam4);
% % nam4 = visualisation3D(y_trial_Tot{2}, 6, totalTimeTrial(2), 0,[0, 0, 1], nam4);
% % nam4 = visualisation3D(y_trial_Tot{3}, 6, totalTimeTrial(3),0,[0, 1, 0], nam4);
% % nam4 = visualisation3D(y_trial{1}, 3, nbData, 0, '.r', nam4);
% % nam4 = visualisation3D(y_trial{2}, 3, nbData, 0, '.b', nam4);
% % nam4 = visualisation3D(y_trial{3}, 3, nbData, 0,'.g', nam4);
% % 
% % title('Trial trajectory');
% % xlabel('x')
% % ylabel('y')
% % zlabel('z')
% % legend('Right','Ahead','Top');
% clear Br Ba Bt i j t yaff