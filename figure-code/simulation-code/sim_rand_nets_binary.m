% sim_rand_nets_binary.m
%
% Simulates the binary network data in Figure 4 of the paper
% "Multistability in neural systems with random cross-connections"
%
% Simulates multiple newtorks across a parameter grid of s and g values.
% This script can be run either locally or on a slurm cluster.

tic


%% Sim options

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select Figure to simulate: fig4
figSim = 'fig4';

% For Figure 4, run with N=10, N=50, and N=100
N = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Determine if HPCC or single PC
if ~isempty(getenv('SLURM_JOB_ID'))
    usingHPCC = 1;
    job_idx = str2double(getenv('SLURM_JOB_ID'));
    task_idx = str2double(getenv('SLURM_ARRAY_TASK_ID'));
    arrayJob_idx = str2double(getenv('SLURM_ARRAY_JOB_ID'));
    slurmNodeList = getenv('SLURM_JOB_NODELIST');
    disp(['slurmNodeList ', slurmNodeList])
    numCores = 1;
    numNets = 1;
else
    usingHPCC = 0;
    numCores = parcluster('local').NumWorkers-1;
    numNets = 100;
end

%% Environment information to record
simulationFileName = mfilename; % Name of this file
simulationFilePath = mfilename('fullpath'); % Path to this file
startDateTime = clock;      % Date and time the script was started
matlabVersion = version;    % version of matlab used
systemVersion = ver;        % Versions of matlab, toolboxes, Java, and OS
[~, git_hash_string] = system('git rev-parse HEAD'); % Git head hash string

%% Simulation options set-up

switch figSim
    case 'fig4'
        useSameNets     = false;    % If true, rng set by ithNet, so identical across (s,g) grid
        runAllTrials    = false;    % If false, only simulates until 2 point attractors are found
        numTrials = 200;            % Number of trials to simulate for each network
        allCombs = 0;               % If true, simulate all 2^N initial conditions (overrides numTrials)
        min_s = 0.0; step_s = 0.2; max_s = 1.4;
        ss = [min_s:step_s:max_s];  % Vector of s values to simulate
        min_g = 0.0; step_g = 0.2; max_g = 3.0;
        gs = [min_g:step_g:max_g];  % vector of g values to simulate
    otherwise
        error('Unkown choice for figSim')
end

if allCombs
    numTrials = 2^N;
end

%% Default parameters
dt = 1;         % Time step (discrete time)
maxT = 20000;   % Max simulation time, to conclude chaos
checkInterval = 50; % How frequently to check for a fixed point or limit cycle
xth = 1.0;      % Threshold for the binary FI curve
pON = 0.5;      % Initial condition probability each unit is ON

%% Names for saving results
if usingHPCC
    dataFolder = strcat('binaryRandDist', ...
        '_N',         num2str(N), ...
        '_nNets',   num2str(numNets), ...
        '_halt',    num2str(~runAllTrials), ...
        '_nTrials', num2str(numTrials), ...
        '_maxT',    num2str(round(maxT/1000)), 'k', ...
        '_sameRNG', num2str(useSameNets), ...
        '_g',         num2str(min_g), '_', ...
        num2str(step_g), '_', ...
        num2str(max_g), ...
        '_s',         num2str(min_s), '_', ...
        num2str(step_s), '_', ...
        num2str(max_s), ...
        '/');
    filename = strcat('binaryRandDist', ...
        '_job', num2str(arrayJob_idx), ...
        '_task', num2str(task_idx,'%05i'));
else
    dataFolder = '../results/';
    filename = strcat('binaryRandDist', ...
        '_N',         num2str(N), ...
        '_nNets',   num2str(numNets), ...
        '_halt',    num2str(~runAllTrials), ...
        '_nTrials', num2str(numTrials), ...
        '_maxT',    num2str(round(maxT/1000)), 'k', ...
        '_sameRNG', num2str(useSameNets), ...
        '_g',         num2str(min_g), '_', ...
        num2str(step_g), '_', ...
        num2str(max_g), ...
        '_s',         num2str(min_s), '_', ...
        num2str(step_s), '_', ...
        num2str(max_s), ...
        '.mat');
end

%% Data matrices
matRNGs     = cell(numNets, numel(ss), numel(gs));
nPointAttrs = nan(numNets, numel(ss), numel(gs));
if runAllTrials
    nTrialsStored = numTrials;
elseif ~runAllTrials
    nTrialsStored = 2;
end
trialFoundMulti = nan(numNets, numel(ss), numel(gs));               % The ith trial on which multistability was found
finalStatesID   = nan(numTrials, numNets, numel(ss), numel(gs));    % Records which state was reached on each trial
integrationTime = nan(numTrials, numNets, numel(ss), numel(gs));    % How long each simulation took to reach a steady state
numUnitsOn      = nan(numTrials, numNets, numel(ss), numel(gs));    % Number of active units in the steady state
stopCriterion   = nan(numTrials, numNets, numel(ss), numel(gs));    % The stopping condition of each trial
% stopCriterion: initialized to nans,
% 1 for t(end)<maxTime, 2 for exp convg. to OFF point attractor,
% 3 for non-convergence within maxT, 4 for binary unit limit cycle

% OFF state is assumed stable for trial=1
stopCriterion(1,:,:,:)   = 1;
numUnitsOn(1,:,:,:)      = 0;
integrationTime(1,:,:,:) = 0;

disp( ['Starting Job: binaryRandDist', ...
    ', N',         num2str(N), ...
    ', NumNets',   num2str(numNets), ...
    ', halt',      num2str(~runAllTrials), ...
    ', NumTrials', num2str(numTrials), ...
    ', allComb',   num2str(allCombs), ...
    ', g',         num2str(min_g), '_', ...
    num2str(step_g), '_', ...
    num2str(max_g), ...
    ', s',         num2str(min_s), '_', ...
    num2str(step_s), '_', ...
    num2str(max_s), newline ])

%% Simulation loop
for si = 1:numel(ss)
    s = ss(si);
    
    for gi = 1:numel(gs)
        g = gs(gi);
        
        disp(['s:', num2str(s, '%.2f'), ' g:', num2str(g, '%.3f'), ...
            ' started at RunTime=', num2str(toc/60/60), ' hours'])
        
        parfor (net = 1:numNets, numCores)
        %for net = 1:numNets
                        
            % Set rng seed for current network
            if usingHPCC
                % Necessary to make sure HPCC workers diverge in the subsequent randnum
                if useSameNets
                    rng(job_idx, 'twister') % All nodes will set same rng
                else
                    rng(mod(job_idx*si*gi, 2^32), 'twister')
                end
            else
                if useSameNets
                    rng(net, 'twister')
                end
            end
            randnum = randi(2^32-1, 1);
            rng(randnum, 'twister')
            rngState = rng;                     % The state of the rng
            matRNGs(net,si,gi) = {rngState};    % Save net's RNG state
            
            % Generate a new random network
            W = sqrt(1/N)*randn(N);
            W = W - diag(diag(W)); W = W';
            W_all = [W*g + s*eye(N)];
            
            
            netFinalStates = nan(nTrialsStored, N); % Accumulate all final states for a given net, then find unique states in it
            netFinalStates(1, :) = zeros(1, N);     % OFF state is always stable
            netStatesFound = 1;
            
            % Simulate network for numTrials
            for trial = 2:numTrials
                
                % Generate initial conditions for this trial
                if allCombs
                    I_init = double([de2bi(trial-1, N, 'left-msb')]' );
                else
                    I_init = ([rand(N,1)<pON]);
                end
                                
                t = 0:dt:checkInterval; % initial time to integrate over
                startT = 2;
                
                x = nan(size(t, 1), N);
                x(1,:) = I_init;
                
                while true 
                    % Continue integrating until either:
                    % 1) a fixed point or limit cycle is found, or
                    % 2) the maxT time is reached
                    
                    for i = startT:numel(t)
                        x(i,:) = [ W_all*x(i-1,:)' ] > xth;
                    end
                    
                    if runAllTrials
                        if x(i,:)==x(i-1,:)
                            % x(i,:) is a fixed point
                            netFinalStates(trial,:) = x(end,:);
                            stopCriterion(trial,net,si,gi) = 1; % ON point attractor
                            attStT = find( all(x(i,:)==x(1:i-1,:), 2));
                            integrationTime(trial, net, si, gi) = attStT(1);
                            numUnitsOn(trial, net, si, gi) = sum(x(i,:));
                            break
                        elseif ~isempty( intersect(x(i,:),x(1:i-1,:),'rows') )
                            % x(i,:) is part of a limit cycle
                            netFinalStates(trial,:) = nan(1, N);
                            integrationTime(trial, net, si, gi) = nan;
                            numUnitsOn(trial, net, si, gi) = nan;
                            stopCriterion(trial,net,si,gi) = 4; % Limit cycle
                            %disp('Limit cycle')
                            break
                        elseif t(end)==maxT
                            netFinalStates(trial,:) = nan(1, N);
                            integrationTime(trial, net, si, gi) = nan;
                            numUnitsOn(trial, net, si, gi) = nan;
                            stopCriterion(trial,net,si,gi) = 3; % Possible chaos
                            %disp('Chaotic?')
                            % x has not reached a fixed point or limit cycle within maxT timesteps
                            break
                        else
                            % continue integrating for another chechInterval number of timesteps
                            t = [t; t(end):dt:(t(end)+checkInterval)];
                            x = [x; nan(checkInterval/dt, N)];
                            startT = startT + checkInterval;
                        end
                        
                    else
                        if x(i,:)==x(i-1,:)
                            % x(i,:) is a fixed point
                            
                            % Trial 1 is OFF attractor, new state found if an ON attractor
                            if sum(x(i,:))>0
                                netStatesFound = netStatesFound+1;
                                netFinalStates(netStatesFound,:) = x(end,:);
                            end
                            if netStatesFound==2
                                trialFoundMulti(net, si, gi) =  trial;
                            end
                            
                            break
                        elseif ~isempty( intersect(x(i,:),x(1:i-1,:),'rows') )
                            %disp('Limit cycle')
                            break
                        elseif t(end)==maxT
                            %disp('Chaotic?')
                            % x has not reached a fixed point or limit cycle within maxT timesteps
                            break
                        else
                            % continue integrating for another chechInterval number of timesteps
                            t = [t; t(end):dt:(t(end)+checkInterval)];
                            x = [x; nan(checkInterval/dt, N)];
                            startT = startT + checkInterval;
                        end
                    end
                    
                    
                end
                
                if netStatesFound==2
                    break
                end
            end % trial loop
            
            % Determine all unique point attractors found across all trials
            if runAllTrials
                [~,~,IC] = unique(netFinalStates, 'rows');
                finalStatesID(:,net, si, gi) = IC;
                nPointAttrs(net, si, gi) = numel(unique(IC(~any(isnan(netFinalStates), 2))));
            else
                nPointAttrs(net, si, gi) = netStatesFound;
            end
            
            
        end     % parfor net loop
        
        %% Save after each paramter point in case of time-limit, overwrites existing file
        if usingHPCC
            [~, ~] = mkdir(dataFolder);
            save([dataFolder, filename])
            % Save run parameters, if first task in sbatch array
            if task_idx==1
                save([dataFolder, 'runParams_job', num2str(arrayJob_idx) ])
            end
        else
            [~, ~] = mkdir(dataFolder);
            save([dataFolder, filename])
        end
        
    end         % g loop
end             % s loop

runTime = toc;


%% Save results
if usingHPCC
    [~, ~] = mkdir(dataFolder);
    save([dataFolder, filename])
    % Save run parameters, if first task in sbatch array
    if task_idx==1
        save([dataFolder, 'runParams_job', num2str(arrayJob_idx) ])
    end
    disp(['Job id: ', num2str(job_idx), newline, ...
        'Task id: ', num2str(task_idx,'%05i'), newline, ...
        'RunTime: ', num2str(runTime/60/60), ' Hours', newline, ...
        'FileName: ', filename, newline])
else
    [~, ~] = mkdir(dataFolder);
    save([dataFolder, filename])
    disp(['RunTime: ', num2str(runTime/60/60), ' Hours', newline, ...
        'FileName: ', filename, newline])
end
