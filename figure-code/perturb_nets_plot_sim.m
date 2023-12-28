% perturb_nets_plot_sim
%
% Code to plot the data produced by perturb_nets_run_sim.m for panel E of 
% Figure 2 of the paper
% "Multistability in neural systems with random cross-connections"

addpath('./functions')

%% Select name and location of data to load

resultsDir = './results/';
simName = '2023-05-28_nNets100_nTrials100_maxT21k_log10pertMag-5_s0_g3.25_FIsig0.1_pert.mat';


%% Load and process data
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle')
load([resultsDir filesep simName])
warning('on','MATLAB:dispatcher:UnresolvedFunctionHandle')

% Analysis options
useFixedBplots = 1;     % Use fits with fixed intercept
rmseType_plots = 'alt'; % Average RMSE over last tauEndAvg of the perturbation
tauEndAvg = 5;          % Taus from end to average over for meanRMSEMat

if useFixedBplots
    fitR2Data = fitR2Mat_fixedb; fitSlopeData = fitSlopetMat_fixedb;
else
    fitR2Data = fitR2Mat; fitSlopeData = fitSlopetMat;
end
switch rmseType_plots
    case 'init'
        RMSEdata = initRMSE;
    case 'end'
        RMSEdata = endRMSE;
    case 'diff'
        RMSEdata = diffRMSE;
    case 'max'
        RMSEdata = maxRMSE;
    otherwise
        RMSEdata = nanmean(meanRMSEMat(:,:,end-round(tauEndAvg*nRMSEBins/linfitend):end), 3);
end

%Parameters for separating chaos from limit cycles
exp_a = 0.025; exp_b = -0.125;
chaos_exp = @(fitR2Data) exp_a*exp(exp_b*fitR2Data);
PtAttRMSDThresh = expectedRMSE/2; % Point attractor if RMSEdata is below this
minPtAttMaxYChange = 1*10^-4;

% Classification: 0=unclassified, 1=OFF Att, 2=ON Att. 3=limit cycle, 4=chaos
scatterColor = zeros(size(fitR2Data));
scatterColor( RMSEdata>chaos_exp(fitR2Data) )    = 4;   % Chaos
scatterColor( RMSEdata<chaos_exp(fitR2Data) )    = 3;   % Limit cycle
scatterColor( RMSEdata<PtAttRMSDThresh & maxUnitChange_y<minPtAttMaxYChange) = 2; % Point attractor
scatterColor( scatterColor==2 & maxUnitVal_y<0.5 ) = 1; % Quiescent point attractor

% Unclassified if chaos or cycle but rates don't vary in original sim (barely stable fixed points?)
scatterColor( [scatterColor==3|scatterColor==4] & maxUnitChange_y<1e-12) = 0;

classLabels = {'Unclas.','Quiescent', 'Stable active', 'Cycle', 'Chaos'};
nClass = numel(classLabels);

%% Net classification, Figure 2E

netDynamicsCount = cellfun(@(c) length(unique(c)), num2cell(scatterColor, 2));
tmp = arrayfun(@(i) ismember([0:(nClass-1)], scatterColor(i,:)), (1:size(scatterColor,1))', 'UniformOutput', false);
netDynamicsBinMat = vertcat(tmp{:});
XX = sortrows(netDynamicsBinMat(:,2:end), 'descend');

% Plot Figure 2E
myPlotSettings(width=2.75, height=1.5)
figure; imagesc(XX)
xticks([1:nClass-1]); xticklabels(classLabels(2:end));
ylabel('Net (sorted)')
colormap(flipud(gray))
ax = gca;
ax.XRuler.Axle.Visible = 'off';
ax.YRuler.Axle.Visible = 'off';
ax.XAxis.TickLength = [0 0];
ax.YAxis.TickLength = [0 0];
myPlotSettings

