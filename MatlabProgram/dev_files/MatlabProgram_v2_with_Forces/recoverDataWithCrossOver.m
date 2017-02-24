%%%%%%%%%%%%%%%%%%%%%%%%%%%% Recover the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% 
% for i=1:nbKindOfTraj
% % %we open the files
% f(i) = fopen(nameD{i}, 'r');
% 
% % we scan the files line per line
% X{i} = textscan(f(i), '%s', 'delimiter', '\n');
% 
% %we close the files
% fclose(f(i));
% 
% %we put data in vectors as numbers
% X{i} = cellfun(@str2num, X{i}{1}, 'UniformOutput', false);
% 
% %we compute the number of block (delimited by a carriage return) and we put
% %data in the vector B.
% numBlock = 1;
% numLine = 0;
% for n = 1:size(X{i}, 1)
%  
%     if isempty(X{i}{n})
%         numBlock = numBlock + 1;
%         numLine = 0;
%     else
%         numLine = numLine+1;
%         data{i}{numBlock}(numLine,:) = X{i}{n}; 
%     end
%  
% end
% 
% end
% clear f numBlock numLine n X ans nameD
 load('Data/dataSere.mat');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create the y vector  and other usefull variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
nbDofTot = size(data{1}{1},2);
%a = 1;
%v = [1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10];
%namf = figure;
%affich=1; %what we want to plot  (x y z fx fy fz))

% % We keep the last sample to try to recognize the movement
for k=1:nbKindOfTraj
    var(k) = size(data{k},2) - 2;
    % totalTime will be the number of time of each trajectories, to begin with it is the same as the totalTime
    % y is the input vector of data
   ymean{k} = zeros(z*(nbDofTot),1);
    for i=1:var(k)
        y{k}{i} = [];
        val = [];
        %to avoid problem after
        totalTime(k,i) =size(data{k}{i},1); 
 
        for j = 1:nbDofTot
           y{k}{i}=  [ y{k}{i} ; data{k}{i}(1:totalTime(k,i),j) ];% data{k}{i}(:,2); data{k}{i}(:,3); filter(v,a,data{k}{i}(:,4)) ; filter(v,a,data{k}{i}(:,5)); filter(v,a,data{k}{i}(:,6))];
           
           val =  [val ; data{k}{i}(1:floor(size(data{k}{i},1)/z):floor(size(data{k}{i},1)/z)*100,j)];
        end
       
         ymean{k}= ymean{k} + (val /var(k));
       % namf = visualisation(y{k}{i}, nbDofTot, totalTime(k,i), affich, [k/nbKindOfTraj, 1 - (k/nbKindOfTraj),0.5*(k/nbKindOfTraj)], namf);hold on; %we draw only the first DOF because too much data  
    end 
    
    ymean{k} = ymean{k}(:,1);
end

% %TRIAL DATA
% %we keep the second last sample of each trajectory to see if we recognize
% and the last for the test of crossover
% %correctly the last movement
for k=1:nbKindOfTraj
 i = size(data{k},2);
 totalTimeTrial2(k)= size(data{k}{i},1); %size(data{k}{i},1); %if we continue the movement to the end
 totalTimeTrial(1,k) = size(data{k}{i-1},1) ; %size(data{k}{i},1); %if we continue the movement to the end
 y_trial_Tot{1}{k} = [];
 y_trial_Tot2{k} = [];
 for j=1:nbDofTot
    y_trial_Tot{1}{k} =  [ y_trial_Tot{1}{k} ; data{k}{i-1}(1: totalTimeTrial(1,k),j)  ] ; %data{k}{i}(:,2); data{k}{i}(:,3); data{k}{i}(:,4) ; data{k}{i}(:,5); data{k}{i}(:,6)];
    y_trial_Tot2{k} =  [ y_trial_Tot2{k} ; data{k}{i}(1: totalTimeTrial2(k),j)  ] ; %data{k}{i}(:,2); data{k}{i}(:,3); data{k}{i}(:,4) ; data{k}{i}(:,5); data{k}{i}(:,6)];
 end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ADD trials 
%%%NEW TRIAL DATA
%load('Data/fakeDataTrial.mat');




clear nam* affich data i j k;