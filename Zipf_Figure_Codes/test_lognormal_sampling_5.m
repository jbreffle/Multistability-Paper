% test_lognormal_sampling_5.m
%
% generate states with a log-normal distribution of sizes (where size means 
% the volume of the basin of attraction) then sub-sample to mimick use of
% random initial conditions that result in being in one basin with
% probability proportional to the size of that basin, and look
% at distribution of the sampled states by number of times they are visited.
%
%   This code is used in Figure 5 of the paper
% "Multiple stable states in neural circuits with random cross-connections"
% 
% Paul Miller, June 5, 2023. pmiller@brandeis.edu
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First set up the distribution
clear

Nstates = 10000000;         % no. of states

mu = 30;                    % log-mean of Gaussian
sig = 6;                    % s.d. of log-Gaussian

% Now generate state sizes as exponential function of a Gaussian
% distribution with mean of mu and s.d. of sig.
% the array "state_size" is sorted so the biggest states have lowest index
% in the array (this helps later with the cumulative sum random selection)
state_size = randn(1,Nstates);
state_size = sort(state_size,'descend');
state_size = exp( sig*state_size + mu);

%% Simply check the distribution of state sizes 
figure(1)
h = histogram(log10(state_size));

figure(2)
plot(h.BinEdges(1:end-1) + 0.5*h.BinWidth,log10(h.Values));

drawnow

%% Now sample from the distribution
Nsamples = 1000000;                             % number of samples
sampled_states = zeros(1,Nsamples);             

% The following variables allow for sampling evenly across the cumulative
% distribution of state sizes so that each state is chosen with a
% probability proportional to its size
total_size = sum(state_size);
cumsum_state_sizes = cumsum(state_size);
rand_samples = total_size*rand(1,Nsamples);

% Now loop through all the random samples to see which state it
% corresponds to by position in cumulative distribution
for i = 1:Nsamples
    r = rand_samples(i);
    j = 1;
    while ( r > cumsum_state_sizes(j) )
        j = j+1;
    end
    sampled_states(i) = j;          % record the state j in that sample, i
end

% now sort the states so identical states are neighboring in the array
sampled_states = sort(sampled_states);  

sampled_n = [];
state_number = 1;
index = 1;
index_old = 0;
while ( index < Nsamples ) 
    state_id = sampled_states(index);   % A new state is found
    % then record how many times that state was visited
    sampled_n(state_number) = length(find(sampled_states == state_id) );
    % update index to the next state visited
    index = index + sampled_n(state_number);
    state_number = state_number + 1;    % next state has a new state number
    
    % the next section indicates progress through the array of visits
    if ( index - index_old > 1000 ) 
        disp([index index/Nsamples])
        index_old = index;
    end
end

% next section removes zeros as it disrupts the log, 
% but there should not be any
sampled_n = sampled_n(find(sampled_n > 0 ) );

% the main variable is "sampled_n" which is simply a list of selected states
% with the number of times that state was selected
% Ordering is by actual size of the state not by the number of visits
save('test_lognormal_samples_5.mat')



