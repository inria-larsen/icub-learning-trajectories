   for t=1:nbData     
        K2 = sigma_new2*PSI_update2{t}' * inv(accuracy*eye(size(PSI_update2{t}*sigma_new2*PSI_update2{t}')) + PSI_update2{t}*sigma_new2*PSI_update2{t}');
        mu_new2 = mu_new2 + K2* (ynew2{t}' - PSI_update2{t}*mu_new2);
        sigma_new2 = sigma_new2 - K2*(PSI_update2{t}*sigma_new2);
        nameFig = visualisation(PSI_z*mu_new, 6, z, i, '+g', nameFig);
    end