%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Recover the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%we open the files
fahead = fopen('Data/dataAvant.txt', 'r');
ftop = fopen('Data/dataHaut.txt', 'r');
fRight = fopen('Data/dataRight.txt', 'r');

% we scan the files line per line
Xa = textscan(fahead, '%s', 'delimiter', '\n');
Xt = textscan(ftop, '%s', 'delimiter', '\n');
Xr = textscan(fRight, '%s', 'delimiter', '\n');

%we close the files
fclose(fahead);
fclose(ftop);
fclose(fRight);
 
%we put data in vectors as numbers
Xa = cellfun(@str2num, Xa{1}, 'UniformOutput', false);
Xt = cellfun(@str2num, Xt{1}, 'UniformOutput', false);
Xr = cellfun(@str2num, Xr{1}, 'UniformOutput', false);

%we compute the number of block (delimited by a carriage return) and we put
%data in the vector B.
numBlock = 1;
numLine = 0;
for n = 1:size(Xa, 1)
 
    if isempty(Xa{n})
        numBlock = numBlock + 1;
        numLine = 0;
    else
        numLine = numLine+1;
        Ba{numBlock}(numLine,:) = Xa{n}; 
    end
 
end

numBlock = 1;
numLine = 0;
for n = 1:size(Xt, 1)
 
    if isempty(Xt{n})
        numBlock = numBlock + 1;
        numLine = 0;
    else
        numLine = numLine+1;
        Bt{numBlock}(numLine,:) = Xt{n}; 
    end
 
end


numBlock = 1;
numLine = 0;
for n = 1:size(Xr, 1)
 
    if isempty(Xr{n})
        numBlock = numBlock + 1;
        numLine = 0;
    else
        numLine = numLine+1;
        Br{numBlock}(numLine,:) = Xr{n}; 
    end
 
end

clear fRight ftop fahead numBlock numLine n Xt Xa Xr ans

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create the y vector  and other usefull variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nbKindOfTraj = 3; % We learn three movements

%number of tests for each kind of trajectories
% We keep the last sample to try to recognize the movement
var(1) = size(Br,2) - 1;
var(2) = size(Ba,2) - 1;
var(3) = size(Bt,2) - 1;

% totalTime will be the number of time of each trajectories, to begin with it is the same as the totalTime
% y is the input vector of data
%namf = figure;
%h = nan(var(1) + var(2) + var(3));
for i=1:var(1)
    totalTime(1,i) = size(Br{i},1);
    y{1}{i} = [ Br{i}(:,1) ; Br{i}(:,2); Br{i}(:,3); abs(Br{i}(:,1)*0.01) ; abs(Br{i}(:,2)*0.02); abs(Br{i}(:,3)*0.03)];
%    namf = visualisation(y{1}{i}(1:totalTime(1,i),1), 1, totalTime(1,i), [1, 0, 0], namf, size(namf,2)); %we draw only the first DOF because too much data
end 
for i=1:var(2)
    totalTime(2,i) = size(Ba{i},1);
    y{2}{i} = [ Ba{i}(:,1) ; Ba{i}(:,2); Ba{i}(:,3); abs(Ba{i}(:,1)*0.01) ; abs(Ba{i}(:,2)*0.02); abs(Ba{i}(:,3)*0.03)];
%    namf = visualisation(y{2}{i}(1:totalTime(2,i),1), 1, totalTime(2,i),  [0,0,1], namf, size(namf,2) );
end
for i=1:var(3)
    totalTime(3,i) = size(Bt{i},1);
    y{3}{i} = [ Bt{i}(:,1) ; Bt{i}(:,2); Bt{i}(:,3); abs(Bt{i}(:,1)*0.01) ; abs(Bt{i}(:,2)*0.02); abs(Bt{i}(:,3)*0.03)];
%    namf = visualisation(y{3}{i}(1:totalTime(3,i),1), 1, totalTime(3,i),  [0,1,0], namf, size(namf,2));
end

% %Information about the previous plot that represent each sample of each
% %trajectory
% title('Trajectories used to learn')
% xlabel('time')
% ylabel('Cartesian position of x for each trajectory')
% legend(namf([2 (var(1)+2) (var(1) + var(2) + 2)]),'Right','Ahead','Top');

%TRIAL DATA
%we keep the last sample of each trajectory to see if we recognize
%correctly the last movement
i = size(Br,2);
totalTimeTrial(1) = size(Br{i},1); %if we continue the movement to the end
y_trial_Tot{1} =  [ Br{i}(:,1) ; Br{i}(:,2); Br{i}(:,3); Br{i}(:,1)*0.01 ; Br{i}(:,2)*0.02; Br{i}(:,3)*0.03]; % the desired trajectory
%the first data of the trajectory from where we would have to recognize the movement
y_trial{1} = [ Br{i}(1:nbData,1) ; Br{i}(1:nbData,2); Br{i}(1:nbData,3)];% Br{i}(1:nbData,1)*0.01 ; Br{i}(1:nbData,2)*0.02; Br{i}(1:nbData,3)*0.03]; 

i = size(Ba,2);
totalTimeTrial(2) = size(Ba{i},1);
y_trial_Tot{2} =  [ Ba{i}(:,1) ; Ba{i}(:,2); Ba{i}(:,3); Ba{i}(:,1)*0.01 ; Ba{i}(:,2)*0.02; Ba{i}(:,3)*0.03];
y_trial{2} = [ Ba{i}(1:nbData,1) ; Ba{i}(1:nbData,2); Ba{i}(1:nbData,3)]; %Ba{i}(1:nbData,1)*0.01 ; Ba{i}(1:nbData,2)*0.02; Ba{i}(1:nbData,3)*0.03];

i = size(Bt,2);
totalTimeTrial(3) = size(Bt{i},1);
y_trial_Tot{3} =  [ Bt{i}(:,1) ; Bt{i}(:,2); Bt{i}(:,3); Bt{i}(:,1)*0.01 ; Bt{i}(:,2)*0.02; Bt{i}(:,3)*0.03];
y_trial{3} = [ Bt{i}(1:nbData,1) ; Bt{i}(1:nbData,2); Bt{i}(1:nbData,3)];% Bt{i}(1:nbData,1)*0.01 ; Bt{i}(1:nbData,2)*0.02; Bt{i}(1:nbData,3)*0.03];

% for j=1:nbKindOfTraj 
%     for t=1:nbData
%          for i=1:nbDof(1) %+ nbDof(2)
%              yaff{j}((nbData*(i-1) + t),1) = y_trial_Tot{j}((totalTimeTrial(j)*(i-1) + t),1); %here we have only the coordinate data
%          end
%     end
% end

% % Plot the trials we will use to verify the learning
% name = figure;
% name = visualisation(y_trial_Tot{1}(1:totalTimeTrial(1),1), 1, totalTimeTrial(1), [1, 0, 0], name, size(name,2));hold on;
% name = visualisation(y_trial_Tot{2}(1:totalTimeTrial(2),1), 1, totalTimeTrial(2), [0, 0, 1], name, size(name,2));
% name = visualisation(y_trial_Tot{3}(1:totalTimeTrial(3),1), 1, totalTimeTrial(3), [0, 1, 0], name, size(name,2));
% name = visualisation(y_trial{1}(1:nbData,1), 1, nbData, '.r', name, size(name,2));
% name = visualisation(y_trial{2}(1:nbData,1), 1, nbData, '.b', name, size(name,2));
% name = visualisation(y_trial{3}(1:nbData,1), 1, nbData, '.g', name, size(name,2));
% title('Trial we keep to try the recognition step (one for each kind of trajectory)');
% xlabel('x position');
% ylabel('iteration');
% legend('Right','Ahead','Top');


clear Br Ba Bt i j t yaff