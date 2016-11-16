%In this function, we create basis function matrix corresponding to the
%number of input information we have and the number of basis function we
%have defined with their bandwith h.

function PSI = computeBasisFunction(z,nbFunctions,alpha, totalTime, h)
    %creating the center of basis function model
    for i = 1 : nbFunctions(1) 
        ci(i) = 0.1*(i-1);%*alpha; 
    end
    for i = 1 : nbFunctions(2) 
        ci2(i) = 0.05*(i-1);%*alpha; 
    end

    for t=1:totalTime %z / alpha
        %creating a basis functions model (time*nbFunctions)
        for i = 1 : nbFunctions(1)
            val = -(alpha*t*0.01-ci(i))*(alpha*t*0.01-ci(i))/(2*h);
            basis(t,i) = exp(val);
        end
 
        sumBI = sum(basis(t,:));
        for i = 1 : nbFunctions(1)
            phi(t,i) = basis(t,i) / sumBI;
        end
    end
    for t =1 : totalTime %z / alpha 
        %creating a basis function model for forces, the same as for movements but with M functions
        for i = 1 : nbFunctions(2)
            val = -(alpha*t*0.01-ci2(i))*(alpha*t*0.01-ci2(i))/(h);
            basis2(t,i) = exp(val);
        end
        sumBI = sum(basis2(t,:));
        for i = 1 : nbFunctions(2)
            phi2(t,i) = basis2(t,i) / sumBI;
        end
    end
    
%      %draw the basis function
%     figure;
%      for i=1:nbFunctions(1)
%          plot(phi(:,i),'r'); hold on;
%     %    plot(basis(:,i),'r'); hold on;
%      end
%     for i=1:nbFunctions(2)
%         plot(phi2(:,i),'b'); hold on;
%        % plot(basis(:,i),'r'); hold on;
%     end
%     title('representation of the basis function used for each type of data (joints in red, forces in blue)')
%     xlabel('time')
%     ylabel('basis normalized')
    
    %TODO ameliorate here to pu as much as we have trajectories!
    %CREATING THE MATRIX BLOCK FOR ALL DOF
    PSI = blkdiag(phi,phi,phi, phi2,phi2,phi2);% we have to put nbDOF time the phi matrix 
end
