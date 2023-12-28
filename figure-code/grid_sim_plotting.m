% grid_sim_plotting.m
%
% Code to plot the data produced by sim_rand_nets.m and 
% sim_rand_nets_binary.m for Figure 4 of the paper
% "Multistability in neural systems with random cross-connections"

addpath('./functions')

%% Select location of data
dataPath = './results/';
addpath(dataPath)

%% Select which data sets to plot:

% tanh
files{101}='genTanh_N10_nNets100_halt1_nTrials200_maxT10k_sig1_sameRNG0_g0_0.2_3_s0_0.2_1.4.mat';
files{102}='genTanh_N50_nNets100_halt1_nTrials200_maxT10k_sig1_sameRNG0_g0_0.2_3_s0_0.2_1.4.mat';
files{103}='genTanh_N100_nNets100_halt1_nTrials200_maxT10k_sig1_sameRNG0_g0_0.2_3_s0_0.2_1.4.mat';

% logistic FIsigma=0.2
files{104}='genLogistic_N10_nNets100_halt1_nTrials200_maxT10k_sig0.2_sameRNG0_g0_0.2_3_s0_0.2_1.4.mat';
files{105}='genLogistic_N50_nNets100_halt1_nTrials200_maxT10k_sig0.2_sameRNG0_g0_0.2_3_s0_0.2_1.4.mat';
files{106}='genLogistic_N100_nNets100_halt1_nTrials200_maxT10k_sig0.2_sameRNG0_g0_0.2_3_s0_0.2_1.4.mat';

% logistic FIsigma=0.1
files{107}='genLogistic_N10_nNets100_nTrials200_maxT10k_sig0.1_sameRNG0_g0_0.2_3_s0_0.2_1.4.mat';
files{108}='genLogistic_N50_nNets100_nTrials200_maxT10k_sig0.1_sameRNG0_g0_0.2_3_s0_0.2_1.4.mat';
files{109}='genLogistic_N100_nNets100_halt1_nTrials200_maxT10k_sig0.1_sameRNG0_g0_0.2_3_s0_0.2_1.4.mat';

% binary
files{110}='binaryRandDist_N10_nNets100_nTrials200_maxT20k_sameRNG0_g0_0.2_3_s0_0.2_1.4.mat';
files{111}='binaryRandDist_N50_nNets100_nTrials200_maxT20k_sameRNG0_g0_0.2_3_s0_0.2_1.4.mat';
files{112}='binaryRandDist_N100_nNets100_nTrials200_maxT20k_sameRNG0_g0_0.2_3_s0_0.2_1.4.mat';

filesToPlot = 101:112;

%% Frac nets multistable across (s,g) grid

tic
for ithFile = filesToPlot
    load(files{ithFile}, 'numTrials', 'stopCriterion', 'numUnitsOn', 'nPointAttrs', 'FIsig', 'ss', 'gs', 'N', 'numNets', 'isTanh')
    
    % Calculate number of point attractor states found for each network
    fracNetsMultistable =  squeeze(mean(nPointAttrs>1, 1));
    
    % Throw warning if any sims on the HPCC failed to complete
    if any(isnan(nPointAttrs), 'all')
        warning('Some nets didnt finish simulating')
        nNan = squeeze(sum(isnan(nPointAttrs), 1));
        N, FIsig, nNan
    end
    
    % Plot fraction of nets that are multistable across the (s, g) grid
    myPlotSettings(width=2, height=1.5)
    figure; imagesc(gs, ss, fracNetsMultistable)
    ax = gca; ax.YDir = 'normal'; ax.YTick = ss; ax.XTick = gs;
    xlabel('g'); ylabel('s')
    caxis([0, 1])
    xticks(0:1:3); yticks(0:0.5:1.5)
    
end
toc

% Plot colorbar for the multistability grid
figure; imagesc(gs, ss, []); caxis([0, 1])
cb = colorbar; clabel = 'Fraction of nets multistable'; ylabel(cb, clabel);

myPlotSettings

