function error = find_sig_logistic3(sig,x_vals,eta_vals,mu,alpha,s,g)
% find_sig_logistic3 calculates the error between the assumed value of the
% variance of inputs and the output value that produces variance of inputs
% assuming a logistic f-I curve in the infinite model.

% given that sig is a s.d. it cannot be negative in a self-consistent
% solution, so this simply "redirects": the solver back toward zero
if ( sig < 0 )
    evar = - sig*sig;
    cv = 0;    
else
    % define range of eta for numerical integration
    eta_min = min(eta_vals);
    eta_max = max(eta_vals);
    
    evar = sig*sig;         % input variance of eta field
    
    if ( sig > 0 )          % true unless sig == 0
        
        d_eta = eta_vals(2)-eta_vals(1);        % integration step-size
        
        % eta has a Gaussian distribution
        p_of_eta = exp(-(eta_vals.*eta_vals)/(2*evar) )/sqrt(2*pi*evar);
        
        % Simple (linear) numerical integration to obtain mean of 
        % firing-rate-squared given the input distribution of eta
        cv1 = d_eta*sum(p_of_eta.*phi_logistic(x_vals,mu,alpha).^2);
        
        % erf_low and erf_high are edge-effect factors to account for input
        % beyond the numerical integration range: both are needed for
        % normalization of the integral to 1, but only erf_high contributes
        % to the variance of firing rates as x(eta) = 1 for high eta but
        % x(eta) = 0 for low eta
        erf_low = 1-0.5*erfc(eta_min/(sqrt(2*evar)));        
        erf_high = 0.5*erfc(eta_max/(sqrt(2*evar)));
        
        pnorm = erf_low + erf_high + d_eta*sum(p_of_eta);
        
        cv = cv1 + erf_high; % don't include erf_low as x(eta) = 0 there
        cv = g*g*cv/pnorm;   % variance of the field has a g-squared factor
        
    else;   % only of sig = 0, so numerics won't work
        
        if ( eta_max < 0 )      % all units are active with zero field input
            cv = g*g;
        else
            if ( eta_min > 0 )  % all units are inactive with zero field input
                cv = 0;
            else;       % all units have rate with zero field input
                index = find(eta_vals(2:end).*eta_vals(1:end-1) < 0 );
                cv = g*g*phi_logistic(x_vals(index),mu,alpha)^2;
            end
        end
    end
end
error = cv - evar;      % error is returned for the solver to look for 0
end

