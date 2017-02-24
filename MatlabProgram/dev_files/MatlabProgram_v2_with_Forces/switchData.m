% here we switch the trial trajectory  with another trajectory used
% previously for the training (to do Leave-one-out cross-validation
% iteration)

clear switchTrajectory timeSwitch alphaSwitch PSI_switch;

switchTrajectory = y_trial_Tot2{i};
timeSwitch = totalTimeTrial2(i);
alphaSwitch = alphaTest2(i);
PSI_switch = PSI_test2{i};

clear y_trial_Tot2{i} totalTimeTrial2(i) alphaTest2(i) PSI_test2{i}

y_trial_Tot2{i} = y{i}{nbCrossOver};
totalTimeTrial2(i) = totalTime(i,nbCrossOver);
alphaTest2(i) = alpha2{i}(nbCrossOver);
PSI_test2{i} = PSI{i}{nbCrossOver};

clear totalTime(i,nbCrossOver) y{i}{nbCrossOver} alpha2{i}(nbCrossOver) PSI{i}{nbCrossOver}


y{i}{nbCrossOver} = switchTrajectory;
totalTime(i,nbCrossOver) = timeSwitch;
alpha2{i}(nbCrossOver) = alphaSwitch;
PSI{i}{nbCrossOver} = PSI_switch;