% same_nets_across_g.m
%
% Code to plot the data produced by sim_rand_nets.m for panels B-F of
% Figure 3 of the paper
% "Multistability in neural systems with random cross-connections"

addpath('./functions')

%% Select location of data
dataPath = './results/';
addpath(dataPath)


%% Set up data filenames and analysis options

% Figure 3 data
files{101} = 'genLogistic_N10_nNets100_halt1_nTrials100_maxT5k_sig0.2_sameRNG1_g1_0.025_25_s0_0.01_0.mat';
files{102} = 'genLogistic_N50_nNets100_halt1_nTrials100_maxT5k_sig0.2_sameRNG1_g1_0.025_25_s0_0.01_0.mat';
files{103} = 'genLogistic_N100_nNets100_halt1_nTrials100_maxT5k_sig0.2_sameRNG1_g1_0.025_25_s0_0.01_0.mat';
files{104} = 'genLogistic_N500_nNets100_halt1_nTrials100_maxT5k_sig0.2_sameRNG1_g1_0.025_25_s0_0.01_0.mat';
files{105} = 'genLogistic_N1000_nNets100_halt1_nTrials100_maxT5k_sig0.2_sameRNG1_g1_0.025_25_s0_0.01_0.mat';

% Check nets in Figure 3 that failed to be multistable at any tested g
% Only 2 N=10 nets failed to reach multistability
files{201} = 'genLogistic_N10_nNets1_halt1_nTrials500_maxT5k_sig0.2_sameRNG1_reRun2604444795_g0.25_0.01_1_s0_0.01_0.mat';
files{202} = 'genLogistic_N10_nNets1_halt1_nTrials200_maxT5k_sig0.2_sameRNG1_reRun2472368426_g0.5_0.01_2_s0_0.01_0.mat';

filesToPlot = [101:105];    % Files to plot for panels 3B-F
checkFailedNetRNGs = false; % Verify that all nets are multistable at some g value
testS = 0; % s value to plot, in case a grid of data is loaded in files{ithFile}
plotOnlyCompletedg = false; % If any net failed to finish simulating for a g value, don't plot that g value
NindsToPlot = [1,3,4]; % Which files{ithFile} to include in panel C

%% Run analysis and plot figures

% Load params from first file for set-up
load(files{filesToPlot(1)}, 'numTrials', 'FIsig', 'ss', 'gs', 'N', 'numNets')

% Initialize data matrices
op_N = nan(numel(filesToPlot), 1);                              % Network size N
op_fracCompletedMultistable = nan(numNets, numel(filesToPlot)); % Fraction of g values each net is multistable for (only considering sucessfully simulated g values)
op_fracNetsMultistable = nan(numel(gs), numel(filesToPlot));    % Fraction of nets multistable at each g value
op_fracNetsMultistable_SEM = nan(numel(gs), numel(filesToPlot));% SEM of op_fracNetsMultistable
op_mean_median_g = zeros(1, numel(filesToPlot));                % Mean over newtorks of each nets median multistable g value
op_mean_median_g_var = zeros(1, numel(filesToPlot));            % Varance of op_mean_median_g
op_g155 = zeros(1, numel(filesToPlot));                         % Fraction of nets multistable at g=1.55
op_g155_nNets = zeros(1, numel(filesToPlot));                   % Count of nets multistable at g=1.55

tic
for ithFile = filesToPlot
    fileInd = find(ithFile == filesToPlot);
    load(files{ithFile}, 'stopCriterion', 'numTrials', 'nPointAttrs', 'FIsig', 'ss', 'gs', 'N', 'numNets')
    
    % Get index of the desired s value testS
    testSind = find( (abs(ss-testS)<0.0001) );
    assert(~isempty(testSind))
    
    % Matrix of networks successfully simulated
    X_completed = squeeze(~isnan(nPointAttrs(:,testSind,:)));
    
    % Probaiblity of multistability across g
    fracNetsMultistable =  squeeze(mean(nPointAttrs>1, 1));
    fracNetsMultistable_SEM = fracNetsMultistable.*numNets.*(1-fracNetsMultistable)./numNets/sqrt(numNets);
    if plotOnlyCompletedg
        fracNetsMultistable(~all(X_completed, 1)) = nan;
    else
        fracNetsMultistable(~any(X_completed, 1)) = nan;
        
    end
    
    % Fraction of nets multistable across g
    op_fracNetsMultistable(:,fileInd) = fracNetsMultistable;
    op_fracNetsMultistable_SEM(:,fileInd) = fracNetsMultistable_SEM;
    
    % Net-wise stablility across g
    X_count = squeeze(nPointAttrs);
    multiStab_bool = squeeze(nPointAttrs)>1;
    tmp = multiStab_bool.*gs; tmp(tmp==0)=nan;
    [~, sortInd] = sort(nansum(tmp, 2)./sum(~isnan(tmp), 2), 'ascend');
    
    % Plot net-wise stability across g, Figure 3B
    myPlotSettings(width=1.5, height=1.5)
    figure; imagesc(gs, 1:numNets, X_count(sortInd,:), 'AlphaData', X_completed(sortInd,:));
    caxis([0, 2]); colormap(parula(3));
    xlabel('g'); ylabel('Net (sorted)')
    yticks([1, round(numNets/2), numNets])
    g_ints = [gs(1):1:gs(end)];
    if numel(g_ints)==25
        xticks([g_ints(1), g_ints(5:5:end)])
    else
        xticks(g_ints)
    end
    axis square
    
    % Mean frac nets multistable at g=1.55
    [~, g155Ind] = min(abs(gs-1.55));
    op_g155([filesToPlot==ithFile]) =  mean(multiStab_bool(:,g155Ind));
    op_g155_nNets([filesToPlot==ithFile]) =  numel(multiStab_bool(:,g155Ind));
    
    % For each net, what is its median multistable g value
    B = multiStab_bool.*gs; B(B==0)=nan;
    netMedianmultistableg = nanmedian(B, 2);
    op_mean_median_g([filesToPlot==ithFile]) = nanmean(netMedianmultistableg);
    op_mean_median_g_var([filesToPlot==ithFile]) = nanvar(netMedianmultistableg);
    
    % Stats over X_count
    X_bool = X_count>1;
    fracCompletedMultistable = sum(X_bool, 2) ./ sum(X_completed, 2);
    
    op_fracCompletedMultistable(:,fileInd) = fracCompletedMultistable;
    op_N(fileInd) = N;
    
    % Check RNG of nets that didn't have any multistability
    disp(['Frac multistable at N=', num2str(N), ' is ', num2str(mean(any(multiStab_bool, 2)))])
    if checkFailedNetRNGs && mean(any(multiStab_bool, 2))<1
        load(files{ithFile}, 'matRNGs')
        failedNetInds = find([any(multiStab_bool, 2)==0]);
        disp('RNG seeds of non-multistable nets')
        for i = 1:numel(failedNetInds)
            disp(['The ', num2str(find(sortInd==failedNetInds(i))), '-th net has rng seed ', num2str(matRNGs{failedNetInds(i),1,1}.Seed)])
        end
    end
    
end
toc

% Plot empty figure with colorbar for panel B
myPlotSettings(width=1.5, height=1.5)
figure;
cb = colorbar; caxis([0, 2]); colormap(parula(3));
clabel = 'Number of Point Att.'; ylabel(cb, clabel); cb.YTick=[0.33,1,1.67]; cb.YTickLabel={'0', '1', '2+'};

% Figure 2C, frac g multistable
myPlotSettings(width=1.5, height=1.5)
figure; swarmchart(repmat(1:numel(op_N(NindsToPlot)), size(op_fracCompletedMultistable, 1), 1), op_fracCompletedMultistable(:,NindsToPlot), ...
    15, ...
    'blue', ...
    'filled', ...
    'MarkerFaceAlpha',0.25, ...
    'MarkerEdgeAlpha',0.25, ...
    'XJitter','density')
xticks([1:numel(op_N(NindsToPlot))]); xticklabels([num2str(op_N(NindsToPlot))])
xlabel('N'); ylabel('Frac. g multistable')
hold on;
y = mean(op_fracCompletedMultistable(:,NindsToPlot), 1);
x = 1:numel(y);
err = std(op_fracCompletedMultistable(:,NindsToPlot), 1)./sqrt(size(op_fracCompletedMultistable(:,NindsToPlot), 1));
errorbar(x, y, err, 'k')
xlim([1-0.6, numel(NindsToPlot)+0.6])

% Figure 2D, frac nets multistable across g values
myPlotSettings(width=2.75, height=1.5);
figure; imagesc(gs, 1:numel(op_N), op_fracNetsMultistable', 'AlphaData', ~isnan(op_fracNetsMultistable'));
yticklabels([op_N])
g_ints = [gs(1):1:gs(end)];
if numel(g_ints)==25
    xticks([g_ints(1), g_ints(5:5:end)])
else
    xticks(g_ints)
end
xlabel('g'); ylabel('N')
cb = colorbar; clabel = 'Frac. multistable'; ylabel(cb, clabel);

% Figure 2E, median g multistable
% Mean over nets of the median g-value at which each net is multistable, across N
myPlotSettings(width=1.5, height=1.5)
figure; errorbar(1:numel(op_N), op_mean_median_g, sqrt(op_mean_median_g_var)/sqrt(numNets))
xticks([1:numel(op_N)]); xticklabels([num2str(op_N)])
xlabel('N'); ylabel({'Median mutistable g', 'Mean\pmSEM'})
yline(1.55, 'k:')
box off

% Figure 2F, Frac multistable at g=1.55 across N
myPlotSettings(width=1.5, height=1.5)
binomialError = (op_g155_nNets.*op_g155.*(1-op_g155))./op_g155_nNets; % binomialVariance / nNets
binomialSEM = sqrt(binomialError)/sqrt(numNets);
figure; errorbar(1:numel(op_N), op_g155, binomialSEM)
xticks([1:numel(op_N)]); xticklabels([num2str(op_N)])
xlabel('N'); ylabel({'Frac. multistable', 'at g=1.55 (\pmSEM)'})
ylim([0, 0.6])
box off

myPlotSettings


%% Check that the two re-ran nets are multistable at the newly tested g values
if checkFailedNetRNGs
    % gs(squeeze(nPointAttrs)>1) % all g values the net is multistable at
    % The 99-th net has rng seed 2604444795, is monostable for g=1:0.025:7.875
    % This net is multistable at gs = [0.74:0.01:0.9], net Ind 22 of files{101}
    load(files{201})
    gs(squeeze(nPointAttrs)>1)
    
    % The 100-th net has rng seed 2472368426, is monostable for g=1:0.025:25
    % This net is multistable at gs = [0.33], net Ind 22 of files{101}
    load(files{202})
    gs(squeeze(nPointAttrs)>1)
end
