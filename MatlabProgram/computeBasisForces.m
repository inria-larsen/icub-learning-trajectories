function PSI = computeBasisCoord(z,nbFunction,alpha, totalTime, h, nbData)
   %creating the center of basis function model
    for i = 1 : nbFunction 
        ci2(i) = 0.05*(i-1); 
    end
  
    for t =1 : z / alpha 
        %creating a basis function model for forces, the same as for movements but with M functions
        for i = 1 : nbFunction
            val = -(alpha*t*0.01-ci2(i))*(alpha*t*0.01-ci2(i))/(h);
            basis2(t,i) = exp(val);
        end
        sumBI = sum(basis2(t,:));
        for i = 1 : nbFunction
            phi2(t,i) = basis2(t,i) / sumBI;
        end
    end

    PSI = blkdiag(phi2(1:nbData,:),phi2(1:nbData,:),phi2(1:nbData,:));% we have to put nbDOF time the phi matrix 
end