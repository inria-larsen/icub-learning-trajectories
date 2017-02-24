
clear all
nbKindOfTraj=1
nameD{1} = 'Data/sit2stand-rigid_init.txt'
for i=1:nbKindOfTraj
% %we open the files
f(i) = fopen(nameD{i}, 'r');

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
        data{i}{numBlock}(numLine,:) = X{i}{n}'; 
    end
 
end

end
clear f numBlock numLine n X ans nameD
data{1,1}{1,1}(:,1) = [] %supprime frame
