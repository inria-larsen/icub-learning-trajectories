function PSI = computeBasisCoord(z,nbFunction,alpha, totalTime, h, nbData)
   %creating the center of basis function model
    for i = 1 : nbFunction 
        ci(i) = 0.1*(i-1);
    end
  
    for t=1:z / alpha
        %creating a basis functions model (time*nbFunctions)
        for i = 1 : nbFunction
            val = -(alpha*t*0.01-ci(i))*(alpha*t*0.01-ci(i))/(2*h);
            basis(t,i) = exp(val);
        end
 
        sumBI = sum(basis(t,:));
        for i = 1 : nbFunction
            phi(t,i) = basis(t,i) / sumBI;
        end
    end

    PSI = blkdiag(phi(1:nbData,:),phi(1:nbData,:),phi(1:nbData,:));% we have to put nbDOF time the phi matrix 
end