% inf_logistic_grid.m
%
% Calculates the possible solutions and their stability in the infinite
% limit of a system of randomly coupled units (zero mean, s.d. g/sqrt(N) )
% such that each unit receives a randomly selected input from a Gaussian
% with s.d. of "g" as well as self-excitation of value "s". 
%
% In this case input-output functions are logistic, and designed such that
% in the absence of external input, any disconnected unit can become
% bistable if s > 1.
%
% This code requires the codes: 
%   findQ_logist.m
%   phi_logistic.m
%   iter_sig_logist2.m
%   find_sig_logistic3.m
% 
% This code was used to generate Figures 6 and 7 in the paper 
% Multistability of neural circuits with random cross-connections.
%
% Paul Miller, June 5, 2023 pmiller@brandeis.edu
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

fix_bifurcation_point = 1;  % fixes bifurcation on y-axis to s = 1;

max_sig_flag = 1;           % uses highest x-values to maximize sigma
test_low_sig = 1;           % also check low sigma solution for stability
Zmax = 6;                   % no. of stds of eta to allow multistability
mu = 0;                     % this will be changed, is "x_th" in the paper
alpha = 0.15;               % inverse steepness of f-I curve, "delta" in the paper
log_g_range = 0;            % can be used to see very large-g results
input = 0;                  % assume no constant external input if zero (just shifts mu)

% filename below is for output 
filename = strcat('inf_logistic_iter_sig_alt15a6_alpha',num2str(alpha),'_input',num2str(input),'.mat')

N_eta_array1 = 5e4;         % Base number of steps for integration over eta      
N_eta_array2 = 5e5;         % Finer grid used when there is a bifurcation

% Next section uses the analytic solution for the logistic function to
% place the single-unit bifurcation at s = 1.
if ( fix_bifurcation_point )
    if ( alpha <= 0.25 )
        mu = 0.5 + sqrt(0.25-alpha) + alpha*log(alpha) ...
            -2*alpha*log(0.5 + sqrt(0.25-alpha));
    end
end

mu = mu - input;        % Any input is like a shift in threshold
disp([alpha, mu])

ss = 0:0.01:1.5;        % Set of values of s for the grid
if ( log_g_range )      % Not used in the paper
    lgg = 0:100;        
    gg = 0.01*exp(log(10)*lgg/10);  % For logarithmic g-scale
else
    gg = [0.01:0.01:3];     % Set of values of g for the grid
end

% Initialize all variable used to store results
% There are 2x4 grids according to the multiple solutions possible for most
% of these values. 
% If there are 2 rows, they correspond to different field solutions with
% different variances of eta.
% The four columns are used in regions where there are bifurcations so
% two values of rate of many units for a given input.
% The columns separately assess solutions with different distributions of
% rates across the continuum to assess whether stable solutions are
% possible.
Qlow_tests = -100*ones(length(ss),length(gg));
Q_tests = -100*ones(length(ss),length(gg),2,4);
eta_std = -ones(length(ss),length(gg),2,4);
high_sig_tests = -ones(length(ss),length(gg),2,4);
multi_eta_tests = -ones(length(ss),length(gg),4);
multi_stable = -ones(length(ss),length(gg),2,4);

% eta_star will contain the bifurcation points in eta closest to zero
eta_star_array = 1e9*Zmax*ones(length(ss),length(gg));
sig_low_array = -ones(length(ss),length(gg));

pick_old_range = 0;             % used if starting from partially complete code
multi_eta_range_min = 0;        % used if starting from partially complete code
multi_eta_range_max = 300;      % used if starting from partially complete code

for i_s = 1:length(ss);         % loop through values of s
    s = ss(i_s);
    if ( s < alpha )            % No bifurcation present
        N_eta_array = N_eta_array1;
    else
        N_eta_array = N_eta_array2; % Bifurcation present means more care needed
    end
    N_x_array = N_eta_array;        % Needed to match x with eta
    
    max_x_dev = 25*alpha;           % Assume logistic function is 0 or 1 this far from x_th
    dx = 2*max_x_dev/N_x_array;     % step-size for xval_array
    xval_array = mu-max_x_dev+dx/2:dx:mu+max_x_dev-dx/2;    % values of x
    % eta_of_x_array has solutions of fixed points in dx/dt equation
    eta_of_x_array = xval_array - s./(1 + exp((mu-xval_array)/alpha) );
    N_xvals = length(xval_array);       % size of x and eta arrays
    
    % limits of the calculated values of eta
    max_eta = eta_of_x_array(end); 
    min_eta = eta_of_x_array(1);
    
    if ( max_eta > min_eta )        % should be true unless alpha = 0
        deta_fine = (max_eta-min_eta)/(N_x_array-1);    % equal steps in eta
        eta_val_array = min_eta:deta_fine:max_eta;      % we need this to get x(eta) from eta(x)
        N_etavals = length(eta_val_array);              % should be same as N_xvals
        
        % The next section will interpolate x(eta) along an equally spaced
        % grid in eta using eta(x) which used an equally spaced grid in x
        % x_of_eta_array has 2 columns to separate out distinct solutions 
        % in the bistable region
        x_of_eta_array = zeros(length(eta_val_array),2);
        x_of_eta_array(1,1) = min(xval_array);
        x_of_eta_array(end,1) = max(xval_array);
        x_of_eta_array(1,2) = min(xval_array);
        x_of_eta_array(end,2) = max(xval_array);
        % First fill out the lower branch by beginning at the lowest values
        % and increasing index in the array 
        ix = 1;
        for i_eta = 2:N_etavals-1
            
            % First find right grid points to use
            while ( ( eta_of_x_array(ix) < eta_val_array(i_eta) ) ...
                    && ( ix < N_xvals ) )
                ix = ix+1;
            end
            
            % Next section does linear interpolation and uses the known
            % value of the logistic function as +1 at large input
            if ( ix < N_xvals )
                
                if ( ix > 1 )
                    frac = ( eta_val_array(i_eta) - eta_of_x_array(ix-1) ) ...
                        /( eta_of_x_array(ix) - eta_of_x_array(ix-1) );
                    
                    x_of_eta_array(i_eta,1) = xval_array(ix-1) + frac*dx;
                else
                    x_of_eta_array(i_eta,1) = eta_val_array(i_eta);
                end
            else;   % large x so add "s"
                x_of_eta_array(i_eta,1) = eta_val_array(i_eta) + s;
            end
            
        end
        
        % Now descend through the indices to get the other branch of the
        % curve in the bistable region
        ix = N_xvals;
        for i_eta = N_etavals-1:-1:2
            
            % First find right grid points to use
            while ( ( eta_of_x_array(ix) > eta_val_array(i_eta) ) ...
                    && ( ix > 1 ) )
                ix = ix-1;
            end
            
            % Next section does linear interpolation and uses the known
            % value of the logistic function as +1 at large input
            if ( ix > 1 )
                if ( ix < N_xvals )
                    frac = (  eta_val_array(i_eta) - eta_of_x_array(ix) ) ...
                        /( eta_of_x_array(ix+1) - eta_of_x_array(ix) );
                    
                    x_of_eta_array(i_eta,2) = xval_array(ix) + frac*dx;
                else; % large x so add "s"
                    x_of_eta_array(i_eta,2) = eta_val_array(i_eta)+s;
                end
            else
                x_of_eta_array(i_eta,2) = eta_val_array(i_eta);
            end
        end
        
    else % if max_eta <= min_eta can only happen if alpha = 0 
        disp('Max eta > min eta')
        max_eta = max(eta_of_x_array);
        min_eta = min(eta_of_x_array);
        deta_fine = (max_eta-min_eta)/(N_x_array-1);
        eta_val_array = min_eta:deta_fine:max_eta;
        N_etavals = length(eta_val_array);
        x_of_eta_array = zeros(length(eta_val_array),2);
        x_of_eta_array(:,1) = eta_val_array;
        x_of_eta_array(:,2) = eta_val_array+s;
    end   % end if max_eta > min_eta
    
    
    % to maximize the variance of network input currents (variance of eta)
    % the absolute values of f(x) should be maximized, which for the
    % logistic function means maximizing x.
    dists = (x_of_eta_array);
    [vals index] = max(dists');
    
    x_of_eta_array_maxsig = nan(1,N_etavals);
    for i_eta = 1:length(eta_val_array)
        x_of_eta_array_maxsig(i_eta) = x_of_eta_array(i_eta,index(i_eta));
    end
    
    % to test the quiescent solution the branch with the minimum x at each
    % value of eta is tested
    [vals index] = min(dists');
    x_of_eta_array_minx = nan(1,N_etavals);
    for i_eta = 1:length(eta_val_array)
        x_of_eta_array_minx(i_eta) = x_of_eta_array(i_eta,index(i_eta));
    end
    
    % to minimize Q the values of x that minimize the gradient of f(x) are
    % key and for the logistic function minimum gradient corresponds to
    % furthest from mu.
    dists = abs(x_of_eta_array - mu);
    [vals index] = max(dists');
    
    x_of_eta_array_minQ = nan(1,N_etavals);
    for i_eta = 1:length(eta_val_array)
        x_of_eta_array_minQ(i_eta) = x_of_eta_array(i_eta,index(i_eta));
    end
    
    
    if ( s >= 4*alpha )     % if there is a bifurcation
        
        % xstar1 and xstar2 are the values of x where d/dx(x-sf(x)) = 0 and
        % eta_star1 and eta_star2 are the corresponding values of eta such that
        % eta = x - s.f(x).
        xstar1 = mu - alpha*log( (s/(2*alpha)) - 1 -sqrt(s/alpha).*sqrt( s/(4*alpha) - 1 ) );
        eta_star1 = xstar1 - 2*alpha./(1-sqrt(1-4*alpha./s));
        xstar2 = mu - alpha*log( (s/(2*alpha)) - 1 +sqrt(s/alpha).*sqrt( s/(4*alpha) - 1 ) );
        eta_star2 = xstar2 - 2*alpha./(1+sqrt(1-4*alpha./s));
        
        % the closest bifurcation point to eta=0 is the most relevant one
        if ( abs(eta_star1) < abs(eta_star2) )
            eta_star = eta_star1;
            xstar = xstar1;
        else
            eta_star = eta_star2;
            xstar = xstar2;
        end
        
    else
        % If there is no bifurcation, set the bifurcation point beyond any
        % region it would be considered for instability or multi-stability
        eta_star = 1e9*Zmax;
        x_star = 1e9*Zmax;
    end
    
    eta_star_array(i_s,:) = eta_star;   % does not depend on g
    
    %% This section is carried out if loading a file and starting from a 
    %  prior result and limits the values of g used for searching for
    %  multiple solutions of eta based on the prior row
    if ( i_s > 1 )
        multi_eta_array = max(squeeze(multi_eta_tests(i_s-1,:,:))' );
        if ( length(find(multi_eta_array > 0 )) > 0 )
            multi_eta_range_min = min(find(multi_eta_array > 0 ));
            multi_eta_range_max = max(find(multi_eta_array > 0 ));
        end
    else
        if ( pick_old_range )
            multi_eta_range_min = 1;
            multi_eta_range_max = 10;
        end
    end
    
    %% Now start looping through values of g
    for i_g = length(gg):-1:1
        
        g = gg(i_g);
        
        % Next little section indicates progress of code
        if ( mod (i_g,1 ) == 0 )
            disp([s,g])
            tic
        end
        
        both_tests = 1;
        Qlow = -1;
        multi_eta_check = 0;
        % the following section loops through the extreme two sets of
        % possible values of x, for given eta, either all at the higher +ve
        % x branch (maxsig) or lower +ve more -ve branch (minQ)
        % it looks at stability of the solution with maximum s.d. of eta
        % (sigval)
        while ( ( both_tests > 0 ) && ( both_tests < 3 ) )
            
            if ( both_tests == 1 )
                x_of_eta_array_full = x_of_eta_array_minQ;
            else
                x_of_eta_array_full = x_of_eta_array_maxsig;
            end
            
            % the function iter_sig_logist2 is used to iterate until one or
            % more self-consistent solutions for the s.d of the field, eta 
            % are found
            [sigval_array bimodal] = iter_sig_logist2(x_of_eta_array_full,eta_val_array,mu,alpha,s,g);
                        
            % Now check if more than one field solution was found
            if ( max(sigval_array) > min(sigval_array) + 1e-6*g )
                many_sigs = 1;
            else
                many_sigs = 0;
            end
            
            % If there are multiple solutions, just consider the first and
            % third as the middle one is unstable
            for sigval_index = 1:2:length(sigval_array)
                
                sigval = sigval_array(sigval_index);
                
                % Record s.d. of eta for that solution
                eta_std(i_s,i_g,(sigval_index+1)/2,both_tests) = sigval;
                
                % Test the stability of the solution
                Q = findQ_logist(sigval,x_of_eta_array_full,eta_val_array,deta_fine,mu,alpha,s,g);
                
                % Test if the s.d. is sufficiently large that units receive
                % input in the bistable region (beyond the bifurcation pt)
                if ( (sigval > abs(eta_star)/Zmax ) )
                    high_sig = 1;
                else
                    high_sig = 0;
                end
                
                % Store those results in the array
                Q_tests(i_s,i_g,(sigval_index+1)/2,both_tests) = Q;
                eta_tests(i_s,i_g,(sigval_index+1)/2,both_tests) = sigval;
                high_sig_tests(i_s,i_g,(sigval_index+1)/2,both_tests) = high_sig;
                
                % Requirement for multiple active stable states with units
                % receiving input between bifurcation points
                if ( ( Q < 1 ) && ( high_sig > 0 ) )
                    multi_stable(i_s,i_g,(sigval_index+1)/2,both_tests) = 1;
                else
                    multi_stable(i_s,i_g,(sigval_index+1)/2,both_tests) = 0;
                end
                
            end
            
            % Record if there are multiple solutions of the field, eta
            if ( many_sigs == 1 )
                multi_eta_tests(i_s,i_g,both_tests) = 1;
            else
                multi_eta_tests(i_s,i_g,both_tests) = 0;
            end
            
            % Update to go to next set of x(eta)
            both_tests = both_tests + 1;
            
            % But exit if there are no bifurcation points
            if ( s < 4*alpha ) % 1-to-1 mapping of eta to x
                both_tests = 0;
            end
            
        end % end loop through two branches of bifurcation region
        
        %% check stability of low activation solution if it is a separate soln
        if ( (test_low_sig == 1 )  ) % Flag to separately test stability of solution with minimum rates
            % Next conditional is zero if there can be only one set of
            % rates and one solution of the field, eta, in which case there
            % is no point looking at separate solutions
            if ( ( s >= 4*alpha) || ( multi_eta_tests(i_s,i_g,1) == 1 ) ) 
                x_of_eta_array_full = x_of_eta_array_minx;      % pick rates of minimum value
                
                % find self-consistent s.d. of field, eta
                [sigval_array bimodal] = iter_sig_logist2(x_of_eta_array_full,eta_val_array,mu,alpha,s,g);
                
                % lowest rates mean lowest s.d. of field
                sigval = min(sigval_array);
                
                sig_low_array(i_s,i_g) = sigval;    % store s.d in array
                
                % may be overkill to look for multiple field solutions
                % again, but just in case
                if ( bimodal )
                    multi_eta_check = 1;
                end
                
                % in case multiple field solutions are already found, record that
                if ( length(sigval_array) > 2 )
                    multi_eta_tests(i_s,i_g,4) = 1;
                end
                
                % loop through all field solutions
                for sigval_index = 1:2:length(sigval_array)
                     
                    sigval = sigval_array(sigval_index);
                    
                    % record s.d. of each field solution
                    eta_std(i_s,i_g,(sigval_index+1)/2,4) = sigval;
                    
                    % check each solution's stability
                    Q = findQ_logist(sigval,x_of_eta_array_full,eta_val_array,deta_fine,mu,alpha,s,g);
                    
                    % test if solution causes units to be beyond a
                    % bifurcation point and hence produce active, potential multistability
                    if ( (sigval > abs(eta_star)/Zmax ) )
                        high_sig = 1;
                    else
                        high_sig = 0;
                    end
                    
                    % record all results in array
                    Q_tests(i_s,i_g,(sigval_index+1)/2,4) = Q;
                    eta_tests(i_s,i_g,(sigval_index+1)/2,4) = sigval;
                    high_sig_tests(i_s,i_g,(sigval_index+1)/2,4) = high_sig;
                    
                    if ( ( Q < 1 ) && ( high_sig > 0 ) ) % if there is active multistability, record it
                        multi_stable(i_s,i_g,(sigval_index+1)/2,4) = 1;
                    else
                        multi_stable(i_s,i_g,(sigval_index+1)/2,4) = 0;
                    end
                    
                end
                
                % Record in a separate array for the low-rate state
                Qlow_tests(i_s,i_g) = Q_tests(i_s,i_g,1,4);
            else;   % If there is only one solution
                Qlow_tests(i_s,i_g) = Q_tests(i_s,i_g,1,1);     % record low-rate state as that solution
            end
        end
        
        %% Check if there is a need to iterate through different var(eta) 
        %  based on the range of multi-stability for previous value of "s" 
        % This is not carried out for all values as it is very
        % time-consuming
        if ( ( i_g >= multi_eta_range_min - 2 ) && ( i_g <= multi_eta_range_max + 1) )
            multi_eta_check = 1;
        end
        
        % If multiple field solutions are already found, no need to check
        % again
        if ( multi_eta_check == 1 )
            if ( max(multi_eta_tests(i_s,i_g,:) == 1 ) )                
               multi_eta_check = 0;
            end
        end
        
        %% Now test if neither the maxsig or minQ solutions are ideal for 
        %  finding multistability based on the distibution of firing rates
        %  when there are bifurcations
        iter_Q = 0;             % try alternative distributions of firing rates to see if stability is possible
        iter_eta = 0;           % try alternative distributions of firing rates to see if multiple field solutions are possible
        if  ( s > 4*alpha )         % if there are bifurcations
            
            % it is possible if in one extreme case the system is
            % stable and in another extreme case there are units beyond the
            % bifurcation point but that solution was not stable
            if ( max(max(multi_stable(i_s,i_g,:,:))) == 0 ) && ...
                    ( min(min(abs( Q_tests(i_s,i_g,:,:)))) < 1 ) && ( max(max( high_sig_tests(i_s,i_g,:,:))) == 1  )
                iter_Q = 1;
            end
            % it is also possible if there are multiple field solutions
            % that the high-variance field may be stable whereas previous
            % distributions produced instability
            if ( max(multi_eta_tests(i_s,i_g,:)) == 1 )
                if ( min(abs(Q_tests(i_s,i_g,2,:))) >= 1 )
                    iter_Q = 1;
                end                
            end
            
            % next section means look for multiple field solutions
            if ( multi_eta_check == 1 )
                iter_eta = 1;
            end
        end
        
        
        %%
        % section to adjust x_of_eta to see if there is a partitioning such
        % that Q<1 and eta_std is significantly > 0 that units receive
        % input beyond a bifurcation point *and* to see if there is a
        % partitioning such that multiple field solutions arise
        % 
        if ( ( iter_Q == 1 ) || ( iter_eta == 1 ) )

            Num_iters = iter_Q + iter_eta;
            iter_start = iter_Q + 2*iter_eta*(1-iter_Q);    % will be 1 unless 2 if only iter_eta is done
            
            % grid of potential solutions for s.d. of field, which will be
            % used as input-field to obtain output-field, with difference
            % being placed in err_vec
            sig0_vec = -g/100:g/100:g;
            Ns0 = length(sig0_vec);
            err_vec = zeros(size(sig0_vec)); % will cross zero when there are solutions
            
            clf
            
            % may need to iterate twice if both a high-sigma low-Q solution
            % is sort as well as a multiple-field solution
            for iter_num = iter_start:iter_start+Num_iters-1
                % the variables with "double" are used to look for field
                % solutions, and relate to partitioning x in different
                % directions from the midpoint if needed
                double_iter_max = 1;
                double_iter = 0;
                doubled_irange = 0;
                doubled_idiff = 0;
                old_zero = 0;
                while( double_iter < double_iter_max )  % only happens once if iter_num == 1
                    double_iter = double_iter + 1;
                    
                    % in each case "xrange" contains the indices for which
                    % there are two different values of x that could be
                    % partitioned differently to help achieve the goal
                    if ( iter_num == 1 )
                        xrange = find(x_of_eta_array_maxsig > x_of_eta_array_minQ + dx);
                    else
                        xrange = find(x_of_eta_array_maxsig > x_of_eta_array_minx + dx);
                        old_bi = 0;
                        any_bi = 0;
                        old_shift = 0;
                    end
                    Nrange = length(xrange);
                    
                    % Begin by partitioning about the midpoint of the range
                    % of alternate values
                    if ( double_iter == 1 )
                        irange = round(Nrange/2);
                        idiff = round(Nrange/2);
                    else;   % If it is the second round, use the suggested alternative partition
                        irange = doubled_irange;
                        idiff = doubled_idiff;
                    end
                    
                    % These three variables used in looking for multiple
                    % field solutions
                    old_irange = 1;
                    change_d = 0;       % becomes 1 if there is a change in direction when halving the partition point
                    loop_num = 0;       % counts number of "halvings" of the range
                    
                    % Iterate by successfully halving the range to find the
                    % splitting point between the maxsig and minQ solutions of
                    % x such that there are multiple solutions of sigma
                    while (idiff > 0)               % becomes zero when no further halving of range is possible
                        loop_num = loop_num + 1;    % number of halvings of range
                        
                        % The values of x will be a combination of the
                        % solutions producing highest rate and either those
                        % yielding minimum Q if stability is sought, 
                        % or those yielding minimum rate if multiple fields
                        % are sought, with the partition point at irange
                        x_of_eta_array_full = x_of_eta_array_maxsig;
                        if ( iter_num == 1 )
                            x_of_eta_array_full(xrange(1:irange)) =  ...
                                x_of_eta_array_minQ(xrange(1:irange));
                        else
                            x_of_eta_array_full(xrange(1:irange)) =  ...
                                x_of_eta_array_minx(xrange(1:irange));
                        end
                        
                        % Given the distribution/partition of x, find the
                        % resulting field given the grid of possible input
                        % field, with err_vec containing the difference
                        for isig0 = 1:Ns0
                            sig0 = sig0_vec(isig0);
                            err_vec(isig0) = find_sig_logistic3(sig0,x_of_eta_array_full,eta_val_array,mu,alpha,s,g);
                        end
                        dsqr_err_vec = err_vec(1:end-1).*err_vec(2:end);    % is <= 0 at a zero crossing
                        zero_cross = find(dsqr_err_vec <= 0);               % indices of zero crossings
                        
                        % Next section used in testing if alternate
                        % partitions may give rise to multiple zero
                        % crossings and hence multiple field solutions
                        if ( loop_num == 1 )
                            old_max = max(err_vec);
                            old_zero_cross = zero_cross(1);
                        end
                        
                        % If there are any zero crossings (there should
                        % always be at least one)
                        if ( length(zero_cross) > 0 )
                            sigval = sig0_vec(max(zero_cross));
                            
                            if ( iter_num == 1 ) % if test is for large sigval and low Q
                                if ( sigval > abs(eta_star)/Zmax )  % if s.d. of field is high enough
                                    old_irange = irange;            % store this partition point
                                    irange = irange + floor(idiff/2);   % see if a lower s.d. is possible
                                else
                                    irange = irange - floor(idiff/2);   % else increase s.d. of field
                                        
                                end
                            else                % if test is for multi_eta
                                % Big section to adjust irange so as to
                                % find if multiple zero crossings are
                                % possible at some value by iteration
                                if ( length(zero_cross) > 2 )   % then there are two stable field solutions
                                    idiff = 0;              % no need to try a new partition
                                    double_iter_max = 0;    % no need to try new iterations
                                    old_irange = irange;    % store this partition point
                                    
                                else;   % this should mean just 1 zero crossing
 
                                    % next few lines check for 
                                    % derivative of the error vector crossing 
                                    % zero and storing those turning points                                    
                                    diff_evec = diff(err_vec);
                                    max_min_vec = diff_evec(1:end-1).*diff_evec(2:end);
                                    max_mins = find(max_min_vec <= 0 );
                                    if ( length(max_mins) > 1 )
                                        bimodal = 1;    % multiple turning points
                                    else
                                        bimodal = 0;    % at most one turning point (not bimodal)
                                    end
                                    
                                    if ( bimodal )      % if there are multiple turning points
                                        % highest one should be a maximum
                                        % and check if it is past the last
                                        % zero crossing
                                        if ( max(max_mins) > zero_cross )
                                            irange = irange - floor(idiff/2);   % increases output sigma
                                            old_shift = -1;     % direction of shift of partition point
                                        else
                                            irange = irange + floor(idiff/2);   % decreases output sigma
                                            old_shift = 1;      % direction of shift of partition point
                                        end
                                        old_bi = 1;             % record that bimodality was found
                                        any_bi = 1;             % record that bimodality was found
                                        
                                    else;                       % only a single max or min
                                        
                                        diff_diff_evec = diff(diff_evec);   % 2nd derivative of error vec
                                        % test for inflexion points where
                                        % 2nd derivative crosses zero
                                        dsqr_diff_diff_evec = diff_diff_evec(1:end-1).*diff_diff_evec(2:end);
                                        inflex = find(dsqr_diff_diff_evec <=0 );    % inflexion point(s)
                                        
                                        % if there is/are inflexion points the
                                        % goal is to place them between the
                                        % zero crossings by changing the
                                        % point of partitioning x
                                        if ( length(inflex) > 1 )
                                            if ( max(inflex) >= max(zero_cross) )
                                                irange = irange - floor(idiff/2);   % increases rates
                                                old_shift = -1;     % direction of change in partition
                                            else
                                               irange =  irange + floor(idiff/2);   % decreases rates
                                               old_shift = 1;       % direction of change in partition
                                            end
                                            
                                        else;           % if no bimodality, no inflexion points
                                            
                                            if  ( any_bi )              % if bimodality has been found
                                                if ( old_bi == 1 )      % and the last partition was bimodal
                                                    % need to shift back in
                                                    % the opposite
                                                    % direction back toward
                                                    % a partition where it
                                                    % was bimodal
                                                    if (old_shift == 1 )
                                                        irange = irange - floor(idiff/2);
                                                        old_shift = -1;
                                                    else
                                                        irange = irange + floor(idiff/2);
                                                        old_shift = 1;
                                                    end
                                                    
                                                else;   % prior partition was not bimodal
                                                    % keep shifting
                                                    % partition in the
                                                    % direction toward
                                                    % where it was once
                                                    % bimodal
                                                    if (old_shift == -1 )
                                                        irange = irange - floor(idiff/2);
                                                        old_shift = -1;
                                                    else
                                                        irange = irange + floor(idiff/2);
                                                        old_shift = 1;
                                                    end
                                                end
                                                
                                            else;   % if no bimodality has been found
                                                
                                                double_iter_max = 2;    % next time start in the other direction of partitioning
                                                
                                                % next lines try to test if
                                                % there is a potential
                                                % switch from a high-sig
                                                % solution to a low-sig
                                                % solution, or vice versa,
                                                % such that between these
                                                % cases there may be a
                                                % partitioning with 2
                                                % solutions and choosing
                                                % the new partition point
                                                % (irange) according to the
                                                % most likely direction
                                                if ( ( max(err_vec) > 4*old_max ) || ( zero_cross(1) > 5*old_zero_cross ) )
                                                    irange = irange + floor(idiff/2);
                                                    old_shift = 1;
                                                    change_d = 1;
                                                else
                                                    if ( (max(err_vec) < 0.25*old_max ) || ( zero_cross(1) < 0.2*old_zero_cross ) )
                                                        irange = irange - floor(idiff/2);
                                                        old_shift = -1;
                                                       change_d =1;
                                                    else
                                                        if ( change_d == 1 )
                                                            irange = irange + old_shift*floor(idiff/2);
                                                        else
                                                            if ( doubled_irange == 0 )
                                                                doubled_irange = irange;
                                                                doubled_idiff = idiff;
                                                            end
                                                            if ( double_iter == 1 )
                                                                irange = irange + floor(idiff/2);
                                                            else
                                                                irange = irange - floor(idiff/2);
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                        old_bi = 0;         % no bimodality found with this partition
                                    end
                                    
                                end
                            end
                        else
                            disp('WARNING zero_cross 0')    % hopefully there is always a zero crossing
                            irange = irange - floor(idiff/2);
                        end
                        
                        % max(err_vec) and lowest zero_cross can change a
                        % lot if there is the possibility of 2 field
                        % solutions
                        old_max = max(err_vec);             % max value of err_vec 
                        old_zero_cross = zero_cross(1);     % lowest val of zero crossing
                        
                        idiff = floor(idiff/2);             % halved to home in on ideal partition
                    end     % move on to reduced idiff
                end         % do second iteration if needed for multiple field solutions
                    
                
                irange = old_irange;            % store the partition
                
                % partition must be within the given range
                if ( irange > Nrange )
                    irange = Nrange;
                end
                if ( irange < 1 )
                    irange = 1;
                end
                
                
                %% Should not be necessary to do this section
                while ( ( length(zero_cross)  == 0 ) && ( irange > 0 ) )
                    disp('WARNING zero_cross 0')
                    
                    irange = irange - 1;
                    x_of_eta_array_full = x_of_eta_array_maxsig;
                    x_of_eta_array_full(xrange(1:irange)) =  x_of_eta_array_minQ(xrange(1:irange));
                    for isig0 = 1:Ns0
                        sig0 = sig0_vec(isig0);
                        err_vec(isig0) = find_sig_logistic3(sig0,x_of_eta_array_full,eta_val_array,mu,alpha,s,g);
                    end
                    dsqr_err_vec = err_vec(1:end-1).*err_vec(2:end);
                    zero_cross = find(dsqr_err_vec <= 0);
                    
                end
                
                %% Proceed now with best set of choices of x to maintain a high-sig solution and look for Q < 1
                x_of_eta_array_full = x_of_eta_array_maxsig;
                if ( iter_num == 1 )    % if looking for a low-Q solution
                    x_of_eta_array_full(xrange(1:irange)) =  ...
                        x_of_eta_array_minQ(xrange(1:irange));
                else;                   % if looking for multiple field solutions
                    x_of_eta_array_full(xrange(1:irange)) =  ...
                        x_of_eta_array_minx(xrange(1:irange));
                end
                
                % Repeat the process of finding the zero crossings first 
                % using a grid of values for s.d. of the field
                for isig0 = 1:Ns0
                    sig0 = sig0_vec(isig0);
                    err_vec(isig0) = find_sig_logistic3(sig0,x_of_eta_array_full,eta_val_array,mu,alpha,s,g);
                end
                dsqr_err_vec = err_vec(1:end-1).*err_vec(2:end);
                zero_cross = find(dsqr_err_vec <= 0);
                
                % sigval_array will store all self-consistent values of the
                % s.d. of the field (at zero crossings)
                sigval_array = zeros(size(zero_cross));
                
                sig_fun = @(sig,x_of_eta_array_full,eta_val_array,mu,alpha,s,g) ...
                    find_sig_logistic3(sig,x_of_eta_array_full,eta_val_array,mu,alpha,s,g);
                
                % then iterate to find the zero between the two values of 
                % sig0_vec that produced a zero crossing
                for isig = 1:length(sigval_array)
                    isig0 = zero_cross(isig);
                    sigmin = sig0_vec(isig0);
                    sigmax = sig0_vec(isig0+1);
                    sigval_array(isig) = ...
                        fzero( @(sig) sig_fun(sig,x_of_eta_array_full,eta_val_array,mu,alpha,s,g), ...
                        [sigmin sigmax ]);
                end
                
                % record if there are multiple zero crossings
                if ( max(sigval_array) > min(sigval_array) + 1e-6*g )
                    many_sigs = 1;
                else
                    many_sigs = 0;
                end
                
                % record if there are multiple zero crossings in the array
                if ( iter_num == 2 )
                    multi_eta_tests(i_s,i_g,both_tests) = many_sigs;
                end
                
                % Now loop through self-consistent solutions and check
                % stability of each one
                for sigval_index = 1:2:length(sigval_array)
                    sigval = sigval_array(sigval_index);
                    
                    % Record field s.d. in array
                    eta_std(i_s,i_g,(sigval_index+1)/2,both_tests) = sigval;
                    
                    % Find the stability coefficient, Q
                    Q = findQ_logist(sigval,x_of_eta_array_full,eta_val_array,deta_fine,mu,alpha,s,g);
                    
                    % If field s.d. is large enough that some units are
                    % beyond a bifurcation point so have multiple rates             
                    if ( (sigval > abs(eta_star)/Zmax ) )
                        high_sig = 1;
                    else
                        high_sig = 0;
                    end
                    
                    % Record all values in arrays
                    Q_tests(i_s,i_g,(sigval_index+1)/2,both_tests) = Q;
                    eta_tests(i_s,i_g,(sigval_index+1)/2,both_tests) = sigval;
                    high_sig_tests(i_s,i_g,(sigval_index+1)/2,both_tests) = high_sig;
                    
                    % Requirement for multiple active stable solutions
                    if ( ( Q < 1 ) && ( high_sig > 0 ) )
                        multi_stable(i_s,i_g,(sigval_index+1)/2,both_tests) = 1;
                    end
                    
                end
                
                % Record if there are multiple field solutions but do not
                % overwrite with zero if they are already found
                if ( many_sigs == 1 )
                    multi_eta_tests(i_s,i_g,both_tests) = 1;
                else
                    if ( multi_eta_tests(i_s,i_g,both_tests) < 1 )
                        multi_eta_tests(i_s,i_g,both_tests) = 0;
                    end
                end
                
            end % end section testing if there is a partition between all solutions on one side or the other of the branch
        end % end loop over number of iters (1 or 2)
        
        %%  Display progress after this value of (s,g) is complete
        if ( mod(i_g,1) == 0 )  % change the "1" to a larger number for less frequent output
            disp([' big_sig lsig: ',num2str(squeeze(high_sig_tests(i_s,i_g,1,:))' ), ...
                '  hsig: ', num2str(squeeze(high_sig_tests(i_s,i_g,2,:))' )])
            disp([' Qval lsig: ',num2str(squeeze(Q_tests(i_s,i_g,1,:))' ),  ...
                ' hsig: ', num2str(squeeze(Q_tests(i_s,i_g,2,:))' )])
            disp([' multi_eta ',num2str(squeeze(multi_eta_tests(i_s,i_g,:))' ), ...
                ' Qlow: ',num2str(Qlow_tests(i_s,i_g) )])
            toc
        end
        
    end
    
    % This saves the data after every set of "g" is complete for a given
    % "s"
    save(filename, 'Q_tests','Qlow_tests','high_sig_tests', ...
        'multi_eta_tests','eta_star_array','sig_low_array','gg','ss',...
        'g','s','alpha','mu','input','Zmax','eta_std','multi_stable')
    
    
end

