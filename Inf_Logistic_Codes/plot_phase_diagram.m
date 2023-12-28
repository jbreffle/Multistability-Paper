% plot_phase_diagram.m
%
% Code to plot the different states as a function of s and g using the
% results from inf_limit_calc.m and its associated function.
%
%   The code uses colormap(turbo) scaled between 1 and 7, with the
%   following key for different types of solution:
%
%   1 = Chaos alone (one unstable solution for the field, eta)
%   1.75 = Quiescence + Chaos (two solutions for the field, eta)
%   2.6 = Quiescence alone (one stable solution for the field, eta)
%   4.5 = Two stable states (two solutions for the field, eta)
%   5.25 = Multiple stable active states (one field solution)
%   6 = Quiescence + multiple stable active states (two field solutions)
%   7 = Chaos + stable active states (two field solutions)
% 
%   This code and variants of it are used to make Figures 6 and 7 in 
%   "Multistability in neural circuits with random cross-connections"
%
%   June 5, 2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

% g_zero_add will add dashed lines on the y-axis if it is otherwise unclear
% that the only solution is quiescent for g = 0 and 0 < s < 1.
g_zero_add = 1; 

% Use appropriate filename output by inf_limit_calc.m
load('inf_logistic_iter_sig_alt15a6_alpha0.15_input0.mat')

[imax jmax] = size(Qlow_tests); % size of grid for the phase diagram

% low_stable_array is a grid with value 1 if the quiescent state is stable
low_stable_array = (sig_low_array*Zmax < eta_star_array).*(Qlow_tests < 1);

% results will hold the numbers (1 to 7) according to the states
results = zeros(imax,jmax);

% Following arrays will indicate which states are present
high_sig_array = zeros(imax,jmax);      % 1 if eta field spreads past a bifurcation
high_stable_array = zeros(imax,jmax);   % 1 if high variance (and rate) state is stable
multi_eta_array = zeros(imax,jmax);     % 1 if there are multiple stable field solutions
bistable_array = zeros(imax,jmax);      % 1 if both field solutions are stable
chaos_array = zeros(imax,jmax);         % 1 if there is an unstable high-Q solution

for i = 1:imax;             % loop through y-axis (s)
    ss_val = exist('ss');
    
    if ( ss_val == 1 )
        s = ss(i);          % s is value of self-excitation
    end
    
    for j = 1:jmax;         % loop through x-axis (g)
        
        eta_star = eta_star_array(i,j);     % bifurcation point closest to zero
        Q_mat = squeeze(Q_tests(i,j,:,:));  % set of Q values for all solutions
        Q_stable = abs(Q_mat) < 1;          % set of stable solutions (2x3)
        
        multi_eta_vec = squeeze(multi_eta_tests(i,j,:));  % multiple field solutions (1x3)  
        sigma_mat = squeeze(eta_std(i,j,:,:));      % s.d. of field, eta (2x3)
        
        if ( alpha > 4*s )                      % if there are no bifurcations
            high_sig = 0;                       % no high rate multiple solutions
        else
            high_sig = (sigma_mat*Zmax > eta_star); % test if field reaches multistable region
        end
        
        low_stable = low_stable_array(i,j);         % if quiescent state is stable
        high_stable_mat = high_sig.*Q_stable;       % if high-rate multiple solutions are stable
        
        unstable_sol = 1;                           % test if any solutions must be unstable
        tests = find( max(Q_mat) >= 0 );            % must have Q > 0 (essentially test all i.c.s)    
        for test = tests
            if ( max(Q_mat(:,test)) >= 1 )          % and that solution cannot shift continuously to a stable one
                unstable_sol = unstable_sol*1;
            else
                unstable_sol = 0;
            end
        end
       
        if ( unstable_sol == 0 )     % test if one branch has only unstable soln
            tests2 = find(multi_eta_vec > 0 );
            if ( length(tests2) > 0 )
                
                for test = tests2'
%                     disp(['test ',num2str(test)])
%                     disp(['next line'])
                    
                    [maxQ index] = max(Q_mat(:,test));
%                     
%                     disp(['max Q ',num2str(maxQ)])
%                     disp(['index ',num2str(index)])
                    
                    unstable_index = (test-1)*2 + (index);
                    if ( max(Q_mat(:,test)) >= 1 )
                        stable_indices = find(abs(Q_mat(:) < 1 ));
                        
                        if ( ( sigma_mat(unstable_index) > max(sigma_mat(stable_indices)) ) ...
                                || ( sigma_mat(unstable_index) < min(sigma_mat(stable_indices))  ) )
                            unstable_sol = 1;
                        end
                        
                    end
                end
            end
            
            % final line to ensure -- if both branches are stable at one
            % point then both solutions are stable
            for test = tests2'
                if ( max(abs(Q_mat(:,test))) < 1 )
                    unstable_sol = 0;
                end
            end
        end
               
        % Now properties of possible solution at one set of parameters is 
        % known, record it in the appropriate arrays
        sigma_array(i,j) = max(max(sigma_mat));
        high_sig_array(i,j) = max(max(high_sig));
        high_stable_array(i,j) = max(max(high_stable_mat));
        multi_eta_array(i,j) = max(multi_eta_vec);
        bistable_array(i,j) = max(sum(Q_stable));
        chaos_array(i,j) = unstable_sol;
        
        %% The following section simply allocates a number in "results" 
        %  according to the properties of the solution for those parameters
        if ( high_stable_array(i,j) == 1 )      % multiple active stable states
            if ( multi_eta_array(i,j) == 1 )    % two field solutions
                if ( unstable_sol == 1 )        % one unstable solution
                    results(i,j) = 7;           % multistable + chaos
                else
                    results(i,j) = 6;           % multistable + quiescent
                end
            else
                results(i,j) = 5.25;            % only multistable active
            end
        else
            if ( low_stable )                   % if quiescent state is stable
                if ( max(multi_eta_vec) < 1 )
                        results(i,j) = 2.6;     % only quiescence
                else
                    if ( unstable_sol == 1 )
                        results(i,j) = 1.75;    % quiescence + chaos
                    else
                        results(i,j) = 4.5;     % quiescence + active stable
                    end
                end
            else
                if ( unstable_sol == 1 )
                    if ( multi_eta_array(i,j) == 1 )
                        disp(' double chaos? ')
                        results(i,j) = 0;       % not seen but perhaps possible
                    else
                        results(i,j) = 1;       % only chaotic solution
                    end
                else
                    disp('Error 2')
                end
            end
            
        end
        
    end
end

%% Now display the results in appropriate color
%
set(0,'DefaultLineLineWidth', 3);
set(0,'DefaultAxesFontSize',20);

figure()


colormap(turbo)
imagesc(results )

set(gca,'YDir','normal')
set(gca,'XTick',[0 100 200 300])
set(gca,'YTick',[0 50 100 150])
set(gca,'XTickLabel',{'0', '1', '2', '3'})
set(gca,'YTickLabel',{'0', '0.5', '1', '1.5'})
xlabel('g')
ylabel('s')

caxis([1 7])
set(gcf,'Units','inches')
set(gcf,'Position',[0 0 5 4])

if ( g_zero_add )       % needed for inf Z limit
    hold on
    map = colormap;
    plot([1 1],[1 100],':','Color',map(round(256*2/7),:),'LineWidth',3 )
end
