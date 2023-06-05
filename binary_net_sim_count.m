% binary_net_sim_count.m
%
% Produce binary networks and test combinations of possible on states for
% stability without simulating but by summing inputs.
%
% Used for the simulations of Figure 9 in the paper 
% Multistability in neural circuits with random cross-connections
% 
% June 5, 2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
tic

ran_seed = 301;         % To allow for reproducible nets
rng(ran_seed);

s = 0.5;                % self-coupling strength
g = 1.2;                % std of cross-coupling (to be normalized by sqrt(N) )

gaussian_conn_flag = 1; % set to 1 for use of randn in connections, otherwise +/- g

Nmax = 12;              % maximum number of units to test
Nvec = [Nmax];          % set of network sizes to consider

num_nets = 10000;       % number of nets at each N
ramp_on = 0;            % set to 1 to gradually ramp down initial constant input
osc_test = 0;           % set to 1 to test for oscillating solutions

frac_multi = zeros(size(Nvec));         % frac of nets multistable at each N
frac_stops = zeros(size(Nvec));         % frac of nets reaching fixed point at each N

mean_mean_Non = zeros(size(Nvec));      % mean # of units on in non-zero fixed points
mean_min_Non = zeros(size(Nvec));       % mean across nets of each N of minimum # units on per net
mean_max_Non = zeros(size(Nvec));       % mean across nets of each N of max # of units on per net
mean_std_Non = zeros(size(Nvec));       % mean across nets of each N of std # units on per net

mean_mean_steps = zeros(size(Nvec));    % mean across nets of each N of mean steps to steady state
mean_max_steps = zeros(size(Nvec));     % mean across nets of each N of max steps to steady state
%% LOOP OVER NETS OF DIFFERENT SIZE
for i_size = 1:length(Nvec);            % loop through sizes of nets
    
    N = Nvec(i_size);                   % N is no. of units
    disp(['N = ', num2str(N) ])
    
    state_array = zeros(N+1,N+1,N+1,2);  % To store details of numbers of states

    low_N = 0;
    if ( N < 21 )                       % for small netss do all starting conditions
        low_N = 1;
    end
    
    % otherwise number of units to start "ON" is drawn uniformly from
    % krange
    krange = [round(Nmax/5):Nmax];
    
    frac_on = zeros(1,num_nets);        % store frac of i.c.s that lead to steady state in each net
    Num_on_states = zeros(1,num_nets);        % store frac of i.c.s that lead to steady state in each net
    mean_Non = zeros(1,num_nets);       % mean number of units ON in non-zero s.s. per net
    min_Non = N*ones(1,num_nets);       % min number of units ON in non-zero s.s. per net
    max_Non = zeros(1,num_nets);        % max number of units ON in non-zero s.s. per net
    std_Non = zeros(1,num_nets);        % std of number of units ON in non-zero s.s.s per net
    
    n_stops = zeros(1,num_nets);        % number of times s.s. reached in each net
    n_osc = zeros(1,num_nets);        % number of times s.s. reached in each net
    mean_steps= zeros(1,num_nets);      % mean number of steps to reach s.s. per net
    max_n_steps = zeros(1,num_nets);    % max number of steps to reach s.s. per net
    
    for i_net = 1:num_nets;             % Loop through different random nets of same size
        if ( mod(i_net,100) == 0 )
            disp([i_net num_nets])
        end
        
        if ( gaussian_conn_flag )
            Wij = ( g/sqrt(N) ) *randn(N);  % Connections from Normal Distribution
        else
            Wij = ( g/sqrt(N) ) * (2*randi(2,N)-3);     % Connections are +/- g/sqrt(N)
        end
        
        Wij = Wij + diag(s*ones(N,1)-diag(Wij));    % set diagonals to s
        Wij0 = Wij + diag(ones(N,1)-diag(Wij));     % initially have diagonals as 1
        
        
        stop_test = 0;                  % will become 1 to stop at steady state
        
        max_steps = N*N;                % maximum # of updates before giving up
        ramp_steps = 0;           % No. of steps to ramp down excitation
        n_record = 10000;            % max period of oscillation to test
        max_starts = 100000;            % maximum no. of initial conditions
        if ( low_N )                    % if starting from all statess
            max_steps = 1;              % probably could be 2
            max_starts = 2^N;           % number of states
        end
        
        Num_states = 0;
        rstate = [];
        starts = 0;                     % no. of initial conditions used
        while ( starts < max_starts )   % in this version try all initial conditions
            
            starts = starts + 1;        % tally number of starts used
            step_num = 0;               % number of updates per simulation
            
            if ( low_N );               % use binary to get pattern if all i.c.s useed
                big_pattern = zeros(1,N);
                pattern = dec2bin(starts-1);
                for i = 1:length(pattern)
                    big_pattern(end+1-i) = eval(pattern(i));
                end
                on_units = find(big_pattern == 1 );
                if ( mod(starts,10000) == 0 )
                    disp ( ['inet = ',num2str(inet),'starts = ',num2str(starts),' % ',num2str(100*starts/max_starts) ] )
                end
                
            else;                       % otherwise randomize the initial condition
                perm_units = randperm(N);
                k = krange(randi(length(krange)));
                on_units = perm_units(1:k); % have k units on
                 if ( mod(starts,1000) == 0 ) || ( N > 5000 )
                     disp ( ['starts = ',num2str(starts),' % ',num2str(100*starts/max_starts), ...
                         ' stops = ',num2str(n_stops(i_net)),' frac_on = ',num2str(frac_on(i_net))] )
                 end
                
            end
            
            k = length(on_units);
            rvec0 = zeros(1,N);
            rvec0(on_units) = 1;
            stable_state = 0;
            stop_test = 0;              % set to 0 initially
            match = 0;
            extra_input = zeros(1,N);
            temp_state = zeros(n_record,N);
            while (stop_test == 0 && step_num < max_steps ) % when to stop each sim
                
                step_num = step_num + 1;                % tally no. of updates
                if ( step_num > 10 )
                    individual_switch = 1;
                else
                    individual_switch = 0;
                end
                
                if ( ramp_on )
                    if ( step_num <= ramp_steps )
                        snow = 1 - (1-s)*step_num/ramp_steps;
                        extra_input(on_units) = snow - s;
                    end
                    
                    inputs = sum(Wij(on_units,:),1) + extra_input;        % total input only from "on" units
                else
                    inputs = sum(Wij(on_units,:),1) ;
                end
                
                rvec = zeros(1,N);                      % rate of all units
                rvec(on_units) = 1;                     % 1 for on units, 0 for others
                % rchange is zero for units stable in their current state
                % given their inputs, otherwise it is multiplied by the
                % difference between the einput and that needed for
                % stability (either OFF or ON)
                rchange = abs(-rvec+heaviside(inputs-1)) .* ...
                    abs(inputs - 1);
                
                [max_change change_unit] = max(rchange); % find unit to switch
                
                if ( individual_switch == 0 )
                    change_unit = find(rchange > 0 );
                end
                
                if ( max_change > 0 )                   % is a unit needs switching
                    rvec(change_unit) = 1-rvec(change_unit);    % 0 to 1 or 1 to 0
                    on_units = find(rvec);              % new set of "on" units
                else
                    stop_test = 1;                      % otherwise all units are stable
                    stable_state = 1;
                end
                
                if ( step_num <= ramp_steps )
                    stop_test = 0;
                end
                
                
                if ( osc_test )
                    if ( step_num > ramp_steps )
                        match = 0;
                        for i = 1:n_record
                            if ( rvec == temp_state(i,:) )
                                match = 1;
                                stop_test = 1;
                            end
                        end
                        ind = mod(step_num-ramp_steps,n_record)+1;
                        temp_state(ind,:) = rvec;
                        
                    end
                end
                
                if ( N >= 1000 )
                    if ( mod(starts,1000) == 0 )
                        if ( ramp_on )
                            if ( ( mod(step_num,100) == 0 ) || ( step_num == 1 ) )
                                disp(['step_num = ',num2str(step_num),' Non = ',num2str(length(on_units)), ...
                                    ' extra_input = ',num2str(snow-s) ] )
                            end
                        end
                    end
                end
                    
            end % next step
            
            %% The next section is book-keeping to accumulate and record
            %  values we may want to look at later
            
            Non = length(on_units);         % # active units in stable states
            
            inputs2 = (inputs - s*rvec0).*(1 + rvec0*(sqrt(k/(k-1))-1) );
            
            if ( k > 1 ) 
            l = length(find(inputs2*sqrt( (k-1)/k) > 1-s ) );
            
            l_m_m = length(find(inputs2 > 1));
            
            state_array(k+1,l-l_m_m+1,l+1,1) = state_array(k+1,l-l_m_m+1,l+1,1)+1;
            if ( stable_state )
                state_array(k+1,l-l_m_m+1,l+1,2) = state_array(k+1,l-l_m_m+1,l+1,2)+1;
            end
            end
            
            if ( stable_state )                % if steady state reached rather than max. no. of steps
                
                if ( ( Non > 0 )  && ( match == 0 ) )            % if the system has active units and is not periodic
                    frac_on(i_net) = frac_on(i_net) + 1;        % add 1 to sims that reach an active state
                    mean_Non(i_net) = mean_Non(i_net) + Non;    % update mean # of active units in states reached
                    if ( Non < min_Non(i_net) )
                        min_Non(i_net) = Non;                   % state with lowest # active for this net
                    end
                    if ( Non > max_Non(i_net) )
                        max_Non(i_net) = Non;                   % state with highest # active for this net
                    end
                    std_Non(i_net) = std_Non(i_net) + Non*Non;  % accumulate to calculate std of # active across states reached
                    
                    % see if it is a new state
                    new_state = 1;
                    istate = 0;
                    while ( istate < Num_states ) && ( new_state == 1 )
                        istate = istate + 1;
                        if rvec == rstate(istate)
                            new_state = 0;
                        end
                    end
                    
                    if ( new_state == 1 )
                        Num_states = Num_states + 1;
                        rstate(Num_states,:) = rvec;
                    end
                    
%                     disp(['i_net = ',num2str(i_net),' starts = ',num2str(starts), ...
%                         ' nsteps = ',num2str(step_num),' Num_on = ',num2str(Non), ' k = ',num2str(k)])
%                     
                end
                n_stops(i_net) = n_stops(i_net) + 1;            % for fraction of stims reaching any s.s.
                
                if ( match == 1 ) 
                    n_osc(i_net) = n_osc(i_net) + 1;
                end
                
                mean_steps(i_net) = mean_steps(i_net) + step_num;   % number of steps to steady state
                if ( step_num  > max_n_steps(i_net) )
                    max_n_steps(i_net) = step_num;              % important to keep "max_steps" above max numbers of steps needed per net
                end
                
            end
            
            if ( N > 200 )
                disp(['Num On = ',num2str(Non), 'step_num',num2str(step_num), ...
                    'stop_test ',num2str(stop_test) ])
            end
        end % next starts means new initial conditions with the same network
        
        disp([i_net frac_on(i_net) Num_states n_osc(i_net)])        % display number of trials that reached steady state for the net
        
        Num_on_states(i_net) = Num_states;
        disp(['frac_on',num2str(frac_on)])

    end % next i_net means new network with same number of units
    
    on_indices = find(frac_on > 0 );        % list of nets with some stable ON states
    
    if ( length(on_indices) > 0 )           % if there are nets with stable ON states look at them
        
        % Note that "frac_on" is actually the total number of initial
        % conditions fo each net that reached a stable ON state
        % we divide by number of initial conditions to obtain a fraction
        % later
        mean_Non(on_indices) = mean_Non(on_indices)./frac_on(on_indices);
        
        std_Non(on_indices) = sqrt( std_Non(on_indices)./frac_on(on_indices) ...
            - mean_Non(on_indices).*mean_Non(on_indices) );
        
        % Now take means over all networks of a given size
        mean_mean_Non(i_size) = mean(mean_Non(on_indices));
        mean_max_Non(i_size) = mean(max_Non(on_indices));
        mean_min_Non(i_size) = mean(min_Non(on_indices));
        mean_std_Non(i_size) = mean(std_Non(on_indices));
        
        % Normalize by number of initial conditions here to get the
        % fraction that result in an ON state
        mean_frac_on(i_size) = mean(frac_on(on_indices))/max_starts;
        mean_num_on_states(i_size) = mean(Num_on_states(on_indices));
        max_num_on_states(i_size) = max(Num_on_states(on_indices));

    end
    
    stop_indices = find(n_stops > 0 );    % nets that reached steady state at least once
    if ( length(stop_indices) > 0 )
        mean_steps(stop_indices) = mean_steps(stop_indices)./n_stops(stop_indices);
        mean_mean_steps(i_size) = mean(mean_steps(stop_indices));
        mean_max_steps(i_size) = mean(max_n_steps(stop_indices));
    end
    
    % frac_mult is fraction of nets that had at least one ON state
    frac_multi(i_size) = length(on_indices)/num_nets;
    
    % frac_stops is fraction of nets that reached any steady state at least
    % once
    frac_stops(i_size) = length(stop_indices)/num_nets;
    
    disp(['frac_multi ', num2str(frac_multi(i_size)), ...
        '   frac_stop ', num2str(frac_stops(i_size)) ]);
    
    disp(['frac_on',num2str(frac_on)])
    
    kon = [2:N];
    plot(kon,sum(sum(squeeze(state_array(3:end,:,:,2)),3),2)/num_nets,'x')


end % next i_size


filename = ['binary_count2_N',num2str(N),'_connflag_',num2str(gaussian_conn_flag), ...
    '_s',num2str(s),'_g',num2str(g),'_rnd',num2str(ran_seed),'.mat'];
save(filename,'-v7.3')


toc
