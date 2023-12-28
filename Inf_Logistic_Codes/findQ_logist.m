function Q = findQ_logist(sigval,x_of_eta_array_full,eta_val_array,deta_fine,mu,alpha,s,g)
%findQ Calculates the value of Q which should be less than 1 if all
%eigenvalues are stable
%
evar = sigval*sigval;           % variance of field, eta

if ( evar >deta_fine ) % p_of_eta is a Gaussian distribution centered on zero
    p_of_eta = exp(-(eta_val_array.*eta_val_array)/(2*evar) )/sqrt(2*pi*evar);
else; % treat as a delta function if variance is small
    p_of_eta = zeros(size(eta_val_array));
    index = find(abs(eta_val_array) <= deta_fine );
    if ( length(index) > 0 )
        p_of_eta(index) = 1/(deta_fine*length(index));
    end
end

% Definition of Q from the paper
Q = g*g*deta_fine*sum(p_of_eta.*(1.0./ ...
    (s - 4*alpha*cosh((x_of_eta_array_full-mu)/(2*alpha)).^2)).^2 );

                % Note that contribution to Q for x-values further from mu are
                % negligible.

end

