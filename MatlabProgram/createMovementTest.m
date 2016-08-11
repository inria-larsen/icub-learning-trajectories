nbData=  floor(z /5);
%phi_begin = phi(1:nbData,:);
%phi2_begin = phi2(1:nbData,:);
%difmoy = zeros(4,1);
 for t=1:nbData
     for i=1:nbDof
         ynew{t}(i,1) = y{5}((z*(i-1) + t),1)- 0.0001;
         ytest(t,i) = y{5}((z*(i-1) + t),1)- 0.0001;
         yaff((nbData*(i-1) + t),1) = y{5}((z*(i-1) + t),1)- 0.0001; % to display the sample
     end
     PSI_test = computeBasisFunction (z,nbFunctions,1, z,h);
     %PSI_begin{t} = blkdiag(phi_begin(t,:), phi_begin(t,:), phi_begin(t,:), phi2_begin(t,:), phi2_begin(t,:), phi2_begin(t,:));
     %calculer le psi reel
end