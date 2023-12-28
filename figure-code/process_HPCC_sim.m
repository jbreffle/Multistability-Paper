% process_HPCC_sim.m
%
% This script takes takes a directory that has simulation results produced
% by the hpcc and combines the data into a single .mat file
%
% Simulations run on the cluster are for a single network at each of the
% specified (s, g) values

%% Select file locations

% Location to save results to
dataFolder = './results/';

% Location of the HPCC results to process
HPCCresultsFolder = './results/genTanh_N50_nNets1_nTrials200_maxT10k_sig1_sameRNG0_g0_0.2_3_s0_0.2_1.4';


%% Set up for processing
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle')
paramFile = dir([targetFile filesep 'runParams_job*.mat']);
load([paramFile.folder, filesep, paramFile.name]);
clear finalStatesID integrationTime numUnitsOn matRNGs stopCriterion tempFinalSTates tempFinalStates IC W allMats
simResultsFilenames = dir([targetFile '/*task*.mat']);
numNets = numel(simResultsFilenames);

%% Initialize final matrices
finalStatesID   = zeros(numTrials, numNets, numel(ss), numel(gs));
integrationTime = zeros(numTrials, numNets, numel(ss), numel(gs));
numUnitsOn      = zeros(numTrials, numNets, numel(ss), numel(gs));
matRNGs         = cell(numNets, numel(ss), numel(gs));
stopCriterion   = zeros(numTrials, numNets, numel(ss), numel(gs));
rmsRates        = zeros(numTrials, numNets, numel(ss), numel(gs));
nPointAttrs     = nan(numNets, numel(ss), numel(gs));
jobRunTime      = nan(numNets);

%% Accumulate data by incrementing through each network's .mat file
tic
for ithNet = 1:numNets
    
    try
        % netDir = dir([dataDir, '*_task', num2str(ithNet), '.mat']);
        netDir = simResultsFilenames(ithNet);
        S = load([netDir.folder, filesep, netDir.name]);
    catch
        netDir = dir([targetFile, '*_task', num2str(ithNet), '_*.mat']);
        S = load([netDir.folder, filesep, netDir.name]);
    end
    
    finalStatesID(:,ithNet,:,:,:) = S.finalStatesID;
    integrationTime(:,ithNet,:,:,:) = S.integrationTime;
    numUnitsOn(:,ithNet,:,:,:) = S.numUnitsOn;
    nPointAttrs(ithNet,:,:,:) = S.nPointAttrs;
    
    if isfield(S, 'rmsRates')
        rmsRates(:,ithNet,:,:,:) = S.rmsRates;
    end
    
    if ~isfield(S, 'matRNGs') % not stored in previous versions of the code
        S.matRNGs.Type = 'threefry';
        S.matRNGs.Seed = job_idx;
        S.matRNGs.State = [];
    end
    matRNGs(ithNet,:,:) = S.matRNGs;
    
    if isfield(S, 'stopCriterion') % used in genRand, not binaryRand
        stopCriterion(:,ithNet,:,:,:) = S.stopCriterion;
    end
    
    if ismember(ithNet, [1:round(numNets/10):numNets])
        disp(['ithNet: ', num2str(ithNet), ' runTime: ', num2str(toc)])
    end
    
    if isfield(S, 'runTime') % Doesn't exist of job ran to time limit
        jobRunTime(ithNet) = S.runTime;
    end
    
end

if ~exist('maxTime') && exist('maxT')
    maxTime = maxT;
end


%% Save in same format as data run on local PC

finalStatesID = single(finalStatesID);
integrationTime = single(integrationTime);
numUnitsOn = single(numUnitsOn);
stopCriterion = single(stopCriterion);
if ~isfield(S, 'stopCriterion'); clear stopCriterion; end
clear S

% Set up filename
tmp=split(netDir.folder, filesep); tmp2=split(tmp{end}, '_'); simType=tmp2{1};
if contains(simType,'binary','IgnoreCase',true)
    filename = strcat(simType, ...
        '_N',         num2str(N), ...
        '_nNets',   num2str(numNets), ...
        '_halt',    num2str(~runAllTrials), ...
        '_nTrials', num2str(numTrials), ...
        '_maxT',    num2str(maxTime/1000), 'k', ...
        '_sameRNG',  num2str(useSameNets), ...
        '_g',         num2str(min_g), '_', ...
        num2str(step_g), '_', ...
        num2str(max_g), ...
        '_s',         num2str(min_s), '_', ...
        num2str(step_s), '_', ...
        num2str(max_s), ...
        '.mat');
else
    filename = strcat(simType, ...
        '_N',         num2str(N), ...
        '_nNets',   num2str(numNets), ...
        '_halt',    num2str(~runAllTrials), ...
        '_nTrials', num2str(numTrials), ...
        '_maxT',    num2str(maxTime/1000), 'k', ...
        '_sig', num2str(FIsig), ...
        '_sameRNG',  num2str(useSameNets), ...
        '_g',         num2str(min_g), '_', ...
        num2str(step_g), '_', ...
        num2str(max_g), ...
        '_s',         num2str(min_s), '_', ...
        num2str(step_s), '_', ...
        num2str(max_s), ...
        '.mat');
end

% Make directory if doesn't exist
[~, ~] = mkdir(dataFolder);

% Save final .mat file
if exist([dataFolder, filename], 'file')
    reset(groot); % myPlotSettings seems to mess up dialog boxes
    promptMessage = sprintf('The file:\n%s\n already exists at \n%s\n Do you want to overwrite it?', ...
        filename, dataFolder);
    titleBarCaption = 'Overwrite?';
    selection = questdlg(promptMessage, titleBarCaption);
    if isequal(selection, 'Yes')
        save([dataFolder, filename], '-v7.3')
    else
        disp('RESULTS NOT SAVED')
    end
else
    save([dataFolder, filename], '-v7.3')
end

runTime = toc;
disp(['Final runTime min ', num2str(runTime/60)])

warning('on','MATLAB:dispatcher:UnresolvedFunctionHandle')
