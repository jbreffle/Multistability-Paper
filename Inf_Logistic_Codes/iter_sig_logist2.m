function [sigval_array bimodal] = iter_sig_logist(x_of_eta_array_full,eta_val_array,mu,alpha,s,g)
% sigval_array should return all self-consistent values of sigma
% bimodal should return 1 if there is more than one max or min of the error
% vector, which would suggest two stable solutions are possible for a
% different distribution of "x"

% sig_fun is needed in the MATLAB zero-crossing finder, fzero, but is 
% identical to "find_sig_logistic3" which returns the difference between
% the input s.d. of the field eta and the s.d. of the resulting field, 
% produced from g^2 multiplied by firing-rate-squared
sig_fun = @(sig,x_of_eta_array_full,eta_val_array,mu,alpha,s,g) ...
    find_sig_logistic3(sig,x_of_eta_array_full,eta_val_array,mu,alpha,s,g);

% Set an array of possible self-constistent values of sigma (must be
% between 0 and g)
sig0_vec = -g/100:g/100:g+g/100;
Ns0 = length(sig0_vec);
err_vec = zeros(size(sig0_vec));

% First, across the range of s.d. of field, eta, find the error between
% input s.d. and output s.d.
for isig0 = 1:Ns0
    sig0 = sig0_vec(isig0);
    err_vec(isig0) = find_sig_logistic3(sig0,x_of_eta_array_full,eta_val_array,mu,alpha,s,g);
end
% find zero crossings of the curve by a non-positive value in the shifted
% product of err_vec
dsqr_err_vec = err_vec(1:end-1).*err_vec(2:end);
zero_cross = find(dsqr_err_vec <= 0);

if ( length(zero_cross) > 0 )   % there should be a zero crossing somewhere
    
    sigval_array = zeros(size(zero_cross));
    
    % For each zero crossing found by the grid search, use the solver to
    % identify the zero precisely within the two neighboring grid points
    for isig = 1:length(sigval_array)
        isig0 = zero_cross(isig);
        sigmin = sig0_vec(isig0);
        sigmax = sig0_vec(isig0+1);
        sigval_array(isig) = ...
            fzero( @(sig) sig_fun(sig,x_of_eta_array_full,eta_val_array,mu,alpha,s,g), ...
            [sigmin sigmax ]);
    end
    
else        % this should not be needed, but ... 
    % otherwise use the solver to find zeros from different starting points
    sig0_array = [-0.1*g 0.01*g 0.1*g 0.5*g g];
    for isig = 1:length(sig0_array)
        sig0 = sig0_array(isig);
        sigval_array(isig) = ...
            fzero( @(sig) sig_fun(sig,x_of_eta_array_full,eta_val_array,mu,alpha,s,g), ...
            sig0);
    end
end

% Section to see if there is more than one max or min
diff_evec = diff(err_vec);      % essentially the derivative of error-vec
% next lines to see if the derivative changes sign
max_min_vec = diff_evec(1:end-1).*diff_evec(2:end);
max_mins = find(max_min_vec <= 0 );
if ( length(max_mins) > 1 )     % if derivative changes sign more than once
    bimodal = 1;            % likely multiple solutions "nearby"
else
    bimodal = 0;            % no multiple solutions possible
end


end

