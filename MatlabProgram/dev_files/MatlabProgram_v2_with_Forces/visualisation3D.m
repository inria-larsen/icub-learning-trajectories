% Function that plot the matrix with the color col1. x is the number of line of the
% matrix, y the number of colonnes, bool=1 if we want forces, 0 if we want
% cartesian position, col1 is the color of the fig, nameFig is the fig

function y = visualisation3D(matrix, x, y, bool, col1, nameFig)

tall= size(nameFig,2);
for i=1:x
    for j=1:y
        val(i,j) = matrix(y*(i-1)+j);
    end
end

if(bool==0)%if we want to plot the coordinate
    if(isa(col1, 'char'))
           nameFig(tall + 1 ) = plot3(val(1,:) ,val(2,:),val(3,:), col1); hold on;
    else
          nameFig(tall+ 1) = plot3(val(1,:) ,val(2,:),val(3,:), 'Color', col1); hold on;
    end
else
    if(isa(col1, 'char'))
           nameFig(tall + 1 ) = plot3(val(4,:) ,val(5,:),val(6,:), col1); hold on;
    else
          nameFig(tall+ 1) = plot3(val(4,:) ,val(5,:),val(6,:), 'Color', col1); hold on;
    end
end


 y = nameFig;
end