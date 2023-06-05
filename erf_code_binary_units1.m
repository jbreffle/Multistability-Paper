% erf_code_binary_units1.m 
%%
%
% This code is used for Figure 9 of "Multistability in neural systems with
% random cross-connections".
%
% This alternative uses the result for the k-th largest value out of N 
% randomly selected Gaussian variables to assess if there are k distinct 
% rows, each with a particular set of k values that sum to greater than the 
% threshold.
%
% OK up to about N=100 but then numerical errors are problematic and vpa is
% needed.
% 
%
% 2023-06-02
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
short_k = 0;        % Given prob is near-zero above k = N/2, this option 
                    % can be set to 1 to speed up code for large N.
figures_on = 0;     % Plot results
part_2 = 0;
no_par = 1;
mu = 0;             % mean of the cross connections

N = 18;             % size of network

g = 1.2;            % s.d. of cross-connections
s = 0.5;            % self-connection strength
        
dx = 0.001;         % step of integration
xvec = [-4:dx:4];   % Range of unit Gaussian for integration
Nx = length(xvec);



kin = [1:N-1];  % all possible number of presynaptic units coactive with a given unit
kon = [2:N];    % the number of units that could be in the ON state together
Nk = length(kon);   
if ( short_k)
    Nk = round(Nk/2);
end

[xmat, kmat] = meshgrid(xvec,kon);              % turn two vectors to a grid

p_of_x = exp(-xmat.*xmat)/sqrt(pi);             % Gaussian in x
p_plus_of_x = erfc(xmat)/2;                     % Prob > x
p_minus_of_x = 1-p_plus_of_x;                   % Prob < x

% Pvals is probability density that kth draw from a Gaussian lies at x.
% So dx*Pvals is probability it lies between x and x+dx.
% Note that "k" in "kth draw" is one greater than row number of Pvals.
Pvals = N*p_of_x.* ...
    [factorial(repmat(N-1, size(kmat)))./ ...
    (factorial(kmat-1).*factorial(repmat(N, size(kmat))-kmat))] .* ...
    (p_plus_of_x.^(kmat-1)).*(p_minus_of_x.^(N-kmat));


disp([s,g])

Xth_on = 1-s+eps;          % minimum sum of input needed to stay on
Xth_off = 1+eps;        % maximum sum of input needed to stay off

theta1 =  (Xth_on-mu)./( g*sqrt(2*kin/N) ); % threshold that a random on-unit is above Xth]
theta2 =  (Xth_off-mu)./( g*sqrt(2*kon/N) ); % threshold that a random off-unit is below Xth]

p_of_k = zeros(Nk+1,N+1);
p_of_k3 = zeros(Nk+1,N+1,N);
Num_of_k3 = zeros(Nk+1,N+1,N);
prob__full = zeros(Nk+1,N+1,N);

Num_of_k = zeros(Nk+1,N+1);
p_of_k_old = zeros(1,Nk);

% Separate calculation for each value of k
for i_k = 1:Nk;                         
    k = kon(i_k)                        % number of active units
    
    if ( theta2(i_k) > theta1(i_k) )     % Should always be the case, otherwise multistability is not possible
        
        % ind1 and ind2 are discrete indices bounding the integral range
        ind1 = min(find(xvec > theta1(i_k) ) );
        ind2 = max(find(xvec < theta2(i_k) ) );
        
        if ( length(ind1) > 0 ) && ( length(ind2) > 0 ) % if range is contained within our grid
            
            r1 = xvec(ind1) - theta1(i_k);  % edge contribution to the integral
            r2 = theta2(i_k) - xvec(ind2);  % edge contribution to the integral
                        
            for m = 0:N;   % m is no. of units with input between theta1 and theta2.
                 
                % l is no. of units with input above theta1 (must be > m and k,
                % cannot be more than N or k+m)
                for l = max(m,k):min(N,k+m);    
                    
                    % comb_factor is number of subsets of N units with given m and l 
                    comb_factor = factorial(N)./(factorial(N-l)*factorial(m)*factorial(l-m));
                    % prob_factor is probability of single subset given the
                    % required N draws from a Gaussian with the numbers
                    % l-m; m; N-l in each range given k active units
                    prob_factor = ( (erfc(theta2(i_k))/2)^(l-m) ) ...
                        *( ( erfc(theta1(i_k))/2 - erfc(theta2(i_k))/2 )^m ) ...
                        *( ( 1- erfc(theta1(i_k))/2 )^(N-l) );
                    
                    % prob_full is the probability of a system existing
                    % with that set of m and l given k units are active
                    prob_full(k+1,m+1, l+1) = comb_factor*prob_factor;
                    
                    % next section indicates number of ways of having k
                    % active units given these subsets of inputs
                    Num_of_k3(k,m+1,l) = nchoosek(m,m-l+k);
                 
                end
                
            end
             
        end
    end
    
end

plot(sum(sum(Num_of_k3.*prob_full(2:end,:,2:end),3),2))






