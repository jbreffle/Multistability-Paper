% sim_rand_nets.m
%
% Simulates the tanh and logistic unit data in Figures 3, 4, and 
% Supplemental Figure 2 of the paper
% "Multistability in neural systems with random cross-connections"
%
% Simulates multiple newtorks across a parameter grid of s and g values.
% This script can be run either locally or on a slurm cluster.

tic

%% Sim options

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select Figure to simulate: fig3, figS2, fig4, reRun
figSim = 'fig3';

% For Figure 3: N=10, N=50, N=100, N=500, and N=1000
% For Supplemental Figure 2: N=10, N=50, N=100, N=500, and N=1000
% For Figure 4: N=10, N=50, and N=100
% For re-running the non-multistable nets in Figure 3: N=10
N = 10; % Network size

% For Figure 3: isTanh=0;FIsig=0.2;
% For Supplemental Figure 2: isTanh=0;FIsig=0.2;
% For Figure 4: isTanh=1;FIsig=1; AND isTanh=0;FIsig=0.1; AND isTanh=0;FIsig=0.1;
% For re-running the non-multistable nets in Figure 3: isTanh=0;FIsig=0.2;
isTanh = false;     % If true, then simulate tanh units, otherwise simulate logistic units
FIsig = 0.2;    % The delta of the FI curve

% For re-running the non-multistable nets in Figure 3: 2604444795 AND 2472368426
% should be 0 otherwise
reRunNetSeed = 0;
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
[~, git_hash_string] = system('git rev-parse HEAD');  % Git head hash string

%% Simulation options set-up

switch figSim
    case {'fig3', 'reRun'}
        useSameNets     = true;     % If true, rng set by ithNet, so identical across (s,g) grid
        runAllTrials    = false;	% If false, only simulates until 2 point attractors are found
        reRunNet        = false;	% Run a single net with a specified RNG seed
        numTrials = 100;            % Number of trials to simulate for each network
        allCombs = 0;               % If true, simulate all 2^N initial conditions (overrides numTrials)
        maxTime = 5 * 1000;         % Max simulation time, to conclude chaos
        min_s = 0.0; step_s = 0.01; max_s = 0.0;
        ss = [min_s:step_s:max_s]; % vector of s values to vary
        min_g = 1.0; step_g = 0.025; max_g = 25.0;
        gs = [min_g:step_g:max_g]; % vector of g values to vary
    case 'figS2'
        useSameNets     = true;     % If true, rng set by ithNet, so identical across (s,g) grid
        runAllTrials    = true;     % If false, only simulates until 2 point attractors are found
        reRunNet        = false;	% Run a single net with a specified RNG seed
        numTrials = 100;            % Number of trials to simulate for each network
        allCombs = 0;               % If true, simulate all 2^N initial conditions (overrides numTrials)
        maxTime = 5 * 1000;         % Max simulation time, to conclude chaos
        min_s = 0.0; step_s = 0.01; max_s = 0.0;
        ss = [min_s:step_s:max_s]; % vector of s values to vary
        min_g = 1.0; step_g = 0.05; max_g = 5.0;
        gs = [min_g:step_g:max_g]; % vector of g values to vary
    case 'fig4'
        useSameNets     = false;	% If true, rng set by ithNet, so identical across (s,g) grid
        runAllTrials    = false;	% If false, only simulates until 2 point attractors are found
        reRunNet        = false;	% Run a single net with a specified RNG seed
        numTrials = 200;            % Number of trials to simulate for each network
        allCombs = 0;               % If true, simulate all 2^N initial conditions (overrides numTrials)
        maxTime = 10 * 1000;        % Max simulation time, to conclude chaos
        min_s = 0.0; step_s = 0.2; max_s = 1.4;
        ss = [min_s:step_s:max_s]; % vector of s values to vary
        min_g = 0.0; step_g = 0.2; max_g = 3.0;
        gs = [min_g:step_g:max_g]; % vector of g values to vary
    otherwise
        error('Unkown choice for figSim')
end

if isequal(figSim, 'reRun')
    numNets=1;
    assert(reRunNetSeed~=0)
else
    assert(reRunNetSeed==0)
end

if allCombs
    if isTanh; numTrials = 3^N;
    else; numTrials = 2^N; end
end

if isTanh; simType = 'genTanh';
else; simType = 'genLogistic'; end

%% Default parameters
RelTol = 1e-6;          % ode45 option, default is 1e-3
AbsTol = RelTol/10;     % ode45 option, default is 1e-6
minDx = RelTol*2;       % Threshold deltaY below which to conclude convergence to point attractor
expConvgMax = RelTol/10; % If all units within transFun(0)+/-expConvgExp at t=tMax, then consider as the zero point attractor
roundFig = 2;           % Number of digits to check for similar states

%% Automatic set up

% Calculate Xth, such that s=1 is the threshold for bistability in an isolated unit
logisticThresh  = 0.5 + sqrt(0.25-FIsig) + FIsig*log(FIsig)-2*FIsig*log(0.5 + sqrt(0.25-FIsig));
tanhThresh =  -1* (sqrt(1-FIsig) -0.5*FIsig*log( (1 + sqrt(1-FIsig))/(1-sqrt(1-FIsig)) ));

logisticFun = @(y) (1./(1+exp( -(y-logisticThresh)/FIsig)));
tanhFun = @(y) (tanh( (y-tanhThresh) / FIsig ) );

invLogisticFun = @(r) logisticThresh - FIsig*log( (1./r) - 1);
invTanhFun = @(r) tanhThresh + FIsig*atanh(r);

if isTanh
    transFun = tanhFun;
    invTransFun = invTanhFun;
    rateBounds = [-1, 1];
    minONVal = abs(0.001); %minimum value for a unit to be considered ON
else
    transFun = logisticFun;
    invTransFun = invLogisticFun;
    rateBounds = [0, 1];
    minONVal = 0.5; %minimum value for a unit to be considered ON
end

% Parameters for random initialization of units' rates
init_mu = 0.5; init_alpha = 0.1;
initSigmoid = @(y) (1./(1+exp( -(y-init_mu)/init_alpha)));


%% Names for saving results
if usingHPCC
    dataFolder = strcat(simType, ...
        '_N',       num2str(N), ...
        '_nNets',   num2str(numNets), ...
        '_halt',    num2str(~runAllTrials), ...
        '_nTrials', num2str(numTrials), ...
        '_maxT',    num2str(round(maxTime/1000)), 'k', ...
        '_sig',     num2str(FIsig), ...
        '_sameRNG', num2str(useSameNets), ...
        '_reRun',   num2str(reRunNetSeed), ...
        '_g',       num2str(min_g), '_', ...
        num2str(step_g), '_', ...
        num2str(max_g), ...
        '_s',         num2str(min_s), '_', ...
        num2str(step_s), '_', ...
        num2str(max_s), ...
        '/');
    filename = strcat('genRandDist_task', ...
        '_job', num2str(arrayJob_idx), ...
        '_task', num2str(task_idx,'%05i'));
else
    dataFolder = '../results/';
    filename = strcat(simType, ...
        '_N',       num2str(N), ...
        '_nNets',   num2str(numNets), ...
        '_halt',    num2str(~runAllTrials), ...
        '_nTrials', num2str(numTrials), ...
        '_maxT',    num2str(maxTime/1000), 'k', ...
        '_sig',     num2str(FIsig), ...
        '_sameRNG', num2str(useSameNets), ...
        '_reRun',   num2str(reRunNetSeed), ...
        '_g',       num2str(min_g), '_', ...
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
rmsRates        = nan(numTrials, numNets, numel(ss), numel(gs));    % The root-mean-squared firing rates at steady state
stopCriterion   = nan(numTrials, numNets, numel(ss), numel(gs));    % The stopping condition of each trial
% stopCriterion: initialized to 0s,
% 1 for t(end)<maxTime, 2 for exp convg. to OFF point attractor,
% 3 for non-convergence within maxT, 4 for binary unit limit cycle

disp( ['Starting Job: ', simType, ...
    ', N',         num2str(N), ...
    ', NumNets',   num2str(numNets), ...
    ', halt',      num2str(~runAllTrials), ...
    ', NumTrials', num2str(numTrials), ...
    ', sig', num2str(FIsig), ...
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
        
        parfor(net = 1:numNets, numCores)
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
            
            if reRunNet
                rng(reRunNetSeed, 'twister')
            else
                randnum = randi(2^32-1, 1);
                rng(randnum, 'twister')
            end
            rngState = rng;                     % The state of the rng
            matRNGs(net,si,gi) = {rngState};    % Save net's RNG state
            
            % Generate a new random network
            W = sqrt(1/N)*randn(N);
            W = W - diag(diag(W));
            W_all = [W*g + s*eye(N)]'; % transpose, to multiply on left side of rate vector
            
            % nTrialsStored numTrials
            netFinalStates = nan(nTrialsStored, N); % accumulate all final states for a given net, then determine unique states
            netStatesFound = 0;
            
             % Simulate network for numTrials
            for trial = 1:numTrials
                
                % Generate initial conditions for this trial
                if allCombs
                    if isTanh
                        r_init = (str2num([dec2base(trial-1, 3, N)]'))-1;
                        %I_init = (double([de2bi(trial-1, N, 'left-msb')]' )*2) - 1;
                    else
                        r_init = double([de2bi(trial-1, N, 'left-msb')]' );
                    end
                    r_init(r_init==rateBounds(1)) = r_init(r_init==rateBounds(1))+1e-4;
                    r_init(r_init==rateBounds(2)) = r_init(r_init==rateBounds(2))-1e-4;
                    I_init = invTransFun( r_init );
                else
                    r_init = (initSigmoid(rand(N,1))-0.5)*diff(rateBounds)+mean(rateBounds);
                    I_init = invTransFun( r_init );
                    if trial==1
                        % Start just a bit away from the zero input state
                        I_init = invTransFun( 1e-4*ones(N,1) );
                    end
                end
                
                % Run simulation
                odefun = @(t, y) ( -y + W_all*transFun(y) );
                xoverFcn = @(t, y) eventfun(y, W_all, minDx, transFun);
                odeOptions = odeset('Events',xoverFcn,'RelTol',RelTol,'AbsTol',AbsTol);
                [t,y] = ode45(odefun, [0, maxTime], I_init, odeOptions);
                
                % check if ode45 converged to a fixed-point
                if runAllTrials
                    % Converged to point attractor criterion for ode45()
                    if t(end)<maxTime
                        netFinalStates(trial,:) = transFun(y(end,:));
                        stopCriterion(trial,net,si,gi) = 1;
                        numUnitsOn(trial,net,si,gi) = sum(abs(transFun(y(end,:)))>minONVal);
                        rmsRates(trial,net,si,gi) = sqrt(mean(transFun(y(end,:)).^2));
                        integrationTime(trial, net, si, gi) = t(end);
                        
                        % Meets criterion for converging to the OFF point attractor
                    elseif all( abs(transFun(y(end,:))-transFun(0)) < expConvgMax )
                        netFinalStates(trial,:) = transFun(y(end,:));
                        stopCriterion(trial,net,si,gi) = 2;
                        numUnitsOn(trial,net,si,gi) = 0;
                        rmsRates(trial,net,si,gi) = sqrt(mean(transFun(y(end,:)).^2));
                        integrationTime(trial, net, si, gi) = t(end);
                        
                        % Fails to meet either point attractor criterion
                        % Therefore, chaotic, limit cycle, or transient chaos
                    else
                        netFinalStates(trial,:) = nan(1, N);
                        stopCriterion(trial,net,si,gi) = 3;
                        numUnitsOn(trial,net,si,gi) = sum(abs(transFun(y(end,:)))>minONVal);
                        rmsRates(trial,net,si,gi) = sqrt(mean(transFun(y(end,:)).^2));
                        integrationTime(trial, net, si, gi) = nan;
                    end
                else
                    % Check if new point attractor
                    if t(end)<maxTime
                        isNewState = ~ismember(round(transFun(y(end,:)), roundFig), round(netFinalStates, roundFig), 'rows');
                        if isNewState
                            netStatesFound = netStatesFound+1;
                            netFinalStates(netStatesFound,:) = transFun(y(end,:));
                        end
                        if netStatesFound==2
                            trialFoundMulti(net, si, gi) =  trial;
                        end
    
                        % Meets criterion for converging to the OFF point attractor
                    elseif all( abs(transFun(y(end,:))-transFun(0)) < expConvgMax )
                        isNewState = ~ismember(round(transFun(y(end,:)), roundFig), round(netFinalStates, roundFig), 'rows');
                        if isNewState
                            netStatesFound = netStatesFound+1;
                            netFinalStates(netStatesFound,:) = transFun(y(end,:));
                        end
                        if netStatesFound==2
                            trialFoundMulti(net, si, gi) =  trial;
                        end

                    else
                        % Fails to meet either point attractor criterion
                        % Therefore, chaotic, limit cycle, or transient chaos
                    end
                end
                
                if netStatesFound==2
                    break
                end
            end % trial loop
            
            % Determine all unique point attractors found across all trials
            if runAllTrials
                [~,~,IC] = unique(round(netFinalStates, roundFig), 'rows');
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


%% Save final results
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


%% Functions
function [x,isterm,dir] = eventfun(y, W_all, minDx, transFun)
% Halts ode45() if the event of x being 0 is detected
% x is the rounded change in y for the next time step
dy = -y + W_all*transFun(y);
x = all( abs(dy) < minDx); %machine precision is (1/eps)~=4.5e+15
isterm = 1;
dir = 0;
end