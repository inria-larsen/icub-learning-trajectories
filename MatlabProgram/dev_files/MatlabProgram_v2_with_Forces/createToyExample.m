% we create simple data to test the toolbox.

%In this part we create fake functions tests to test the program. This part
%will not be used when we will have real data
nbKindOfTraj = 1; % number of kind of trajectories (will be three after)

var = 30*ones(nbKindOfTraj,1); %number of tests for each kind of trajectories

figure
moy =10;
sigma = 7;
bruit = moy + sigma*randn(1,30);

totalTime = ones(1,var(1))*100 + bruit; %will be the number of time of each trajectories, to begin with it is the same as the totalTime

bruitf = moy + 2*rand(1,30);
bruiti = moy + 2*rand(1,30);
% for i =1:12
%     bruiti(i) = bruiti(i)*(power(-1,i));
% end

xf = ones(1, var(1))*10 + bruitf ; %final positions of trajectories analysed
xi = zeros(1, var(1)) + bruiti;

for i=1:var(1)
   val = zeros(totalTime(1,i),1);
    for t=1:totalTime(1,i)
        tovd = (t / totalTime(i));
        val(t) = xi(i) + (xf(i) - xi(i))*(10*power(tovd,3) - 15*power(tovd,4) + 6*power(tovd,5)) -1
    end
    yTest{1}{i} = val;
    plot(yTest{1}{i});hold on;
end
