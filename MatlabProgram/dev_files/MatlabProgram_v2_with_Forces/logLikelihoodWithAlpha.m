%function that compute the log likelihood
% If A positif symetric
% R = chol(A) where R'*R = A
%TODO: add possibiliy to add nbDof (not -0.5 for 2pi)
function log_p = logLikelihoodWithAlpha(x,mu,S, alpha, mu_alpha, sigma_alpha)
   logAlpha = -0.5*log(2*pi*sigma_alpha) - ((alpha - mu_alpha)*(alpha - mu_alpha))/(2*sigma_alpha);

    [Sigma, p] = chol(S);
   if(p ~=0) error('Error in the cholesky decomposition');
   end
    logdetSigma = sum(log(diag(Sigma))); % logdetSigma
    log_p = -0.5*log(2*pi) -0.5*logdetSigma -0.5*(x-mu)'*(S\(x-mu))  + logAlpha;  

end
