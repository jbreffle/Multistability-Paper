% rmsrates_across_g.m
%
% This script creates Figure 11

%% Select which panel to create

% Options to create each panel: a1, a2, a3, b1, b2, b3
figOpt = 'b3'; 
% Panel a: delta=0.2, s=0
% Panel b: delta=0.2, s=0.5
% Column 1: All fixed point trials
% Column 2: All quiescent trials
% Column 3: All stable active trials


%% Set up

%% Select location of data
dataPath = './results/';
addpath(dataPath)

% Plot Settings
myPlotSettings(width=2.0, height=1.5, ttlfsz=1.0)

% Panel a: alpha=0.2, s=0.0
files{101} = 'genLogistic_N10_nNets100_halt0_nTrials25_maxT5k_sig0.2_sameRNG0_g1_0.05_5_s0_0.01_0.mat';
files{102} = 'genLogistic_N100_nNets100_halt0_nTrials25_maxT5k_sig0.2_sameRNG0_g1_0.05_5_s0_0.01_0.mat';
files{103} = 'genLogistic_N1000_nNets100_halt0_nTrials25_maxT5k_sig0.2_sameRNG0_g1_0.05_5_s0_0.01_0.mat';

% Panel b: alpha=0.2, s=0.0
files{201} = 'genLogistic_N10_nNets1000_halt0_nTrials25_maxT5k_sig0.2_sameRNG0_reRun0_g0_0.2_3_s0.5_0.01_0.5.mat';
files{202} = 'genLogistic_N100_nNets1000_halt0_nTrials25_maxT5k_sig0.2_sameRNG0_reRun0_g0_0.2_3_s0.5_0.01_0.5.mat';
files{203} = 'genLogistic_N1000_nNets1000_halt0_nTrials25_maxT5k_sig0.2_sameRNG0_g0_0.2_3_s0.5_0.01_0.5.mat';

if isequal(figOpt(2), '1')
    splitOpt = 5;
elseif isequal(figOpt(2), '2')
    splitOpt = 4;
elseif isequal(figOpt(2), '3')
    splitOpt = 3;
else
    error('Unkonwn figOpt')
end
if isequal(figOpt(1), 'a')
    filesToPlot = [101, 102, 103];
    target_s = 0;
elseif isequal(figOpt(1), 'b')
    filesToPlot = [201, 202, 203];
    target_s = 0.5;
else
    error('Unkonwn figOpt')
end

% Create the dummy legend figure
figure(1234) 

% Plot the inf-N expected curve
plot_mean_r

% Create the colormap
cm = colormap(parula(numel(filesToPlot)));

%% Loop over files to plot
tic
legend_N = {};
for ithFile = filesToPlot
    
    load(files{ithFile}, 'numTrials', 'numUnitsOn', 'stopCriterion', 'rmsRates', 'FIsig', 'ss', 'gs', 'N', 'numNets')
    
    try
        tmp = split( files{ithFile}, '_');
        switch  tmp{1}; case 'genLogistic'; simType = 'log'; case 'genTanh'; simType = 'tanh'; end
    end
    
    % rmsRates, mean over all trials, including non-converging
    meanRMSRates = reshape(nanmean(rmsRates, [1, 2]), numel(ss), numel(gs) );
    % figure; histogram(rmsRates)
    
    % rmsRates, mean over all converging trials
    % rmsRates(stopCriterion==3)=nan;
    % meanRMSRates = reshape(nanmean(rmsRates, [1, 2]), numel(ss), numel(gs) );
    
    % stopCriterion: initialized to 0s,
    % 1 for t(end)<maxTime, 2 for exp convg. to OFF point attractor,
    % 3 for non-convergence within maxT, 4 for binary unit limit cycle
    
    sInd = find( abs(target_s-ss)<0.0001)
    
    
    valid_nonconvg = double([stopCriterion==3]); valid_nonconvg(valid_nonconvg==0)= nan;
    valid_convg = double([stopCriterion~=3]); valid_convg(valid_convg==0)= nan;
    valid_convg_active = valid_convg .* [numUnitsOn>=1]; valid_convg_active(valid_convg_active==0)=nan;
    valid_convg_quiesc = valid_convg .* [numUnitsOn<1]; valid_convg_quiesc(valid_convg_quiesc==0)=nan;
    
    X_nonvonvg =     rmsRates(:,:,sInd,:).*valid_nonconvg(:,:,sInd,:);
    X_activestable = rmsRates(:,:,sInd,:).*valid_convg_active(:,:,sInd,:);
    X_quiescstable = rmsRates(:,:,sInd,:).*valid_convg_quiesc(:,:,sInd,:);
    X_stable = rmsRates(:,:,sInd,:).*valid_convg(:,:,sInd,:);
    %figure; plot(gs, squeeze(nanmean(X_nonvonvg, [1,2,3])))
    %figure; plot(gs, squeeze(nanmean(X_activestable, [1,2,3])))
    %figure; plot(gs, squeeze(nanmean(X_quiescstable, [1,2,3])))
    
    
    if splitOpt==1
        X = rmsRates(:,:,sInd,:);
    elseif splitOpt==2
        X = X_nonvonvg;
    elseif splitOpt==3
        X = X_activestable;
    elseif splitOpt==4
        X = X_quiescstable;
    elseif splitOpt==5
        X = X_stable;
    else
        error()
    end
    
    
    yvals = squeeze(mean(X, [1,2,3], 'omitnan'));
    valid_yval = squeeze(sum(~isnan(X), [1,2,3]))'>0;
    plot(gs(valid_yval), yvals(valid_yval), 'Color', cm(find(ithFile==filesToPlot),:))
    
    disp([simType , ', FIsig', num2str(FIsig), ', N', num2str(N), ', nNets', num2str(numNets), ', nTrials', num2str(numTrials)])
    xlabel('g'); ylabel('RMS rates')
    legend_N{end+1} = num2str(N);
    
    
    ylim([0, 0.675])
    if all(filesToPlot>=100) && all(filesToPlot<200)
        xlim([1, 3])
    end
    
    % Plots how many trials underly the means plotted in the main figure
    % figure; plot(gs, squeeze(sum(~isnan(X), [1,2,3]))); xlabel('g'); ylabel('Count (trials)') 
    
end
title(title_list(splitOpt), 'fontweight', 'normal')
toc



%% Dummy figure with corresponding legend
figure(1234); hold on
plot(0, 0, 'k')
for i = 1:size(cm, 1)
    plot(0,0, 'Color', cm(i,:))
end
leg1 = legend(['\infty', legend_N], 'location', 'best', 'box', 'off');
leg1.ItemTokenSize = [12, 12];

% Pull relevant figure to top
% h = findobj('Type','figure');
% figure(h(end))

