% same_nets_across_g_classes.m
%
% Code to plot the data produced by sim_rand_nets.m for supplementary
% Figure 2 of the paper
% "Multistability in neural systems with random cross-connections"

addpath('./functions')

%% Select location of data
dataPath = './results/';
addpath(dataPath)

%% Select which data sets to plot
files{101} = 'genLogistic_N10_nNets100_nTrials100_maxT5k_sig0.2_sameRNG1_g1_0.05_5_s0_0.01_0.mat';
files{102} = 'genLogistic_N50_nNets100_nTrials100_maxT5k_sig0.2_sameRNG1_g1_0.05_5_s0_0.01_0.mat';
files{103} = 'genLogistic_N100_nNets100_nTrials100_maxT5k_sig0.2_sameRNG1_g1_0.05_5_s0_0.01_0.mat';
files{104} = 'genLogistic_N500_nNets100_nTrials100_maxT5k_sig0.2_sameRNG1_g1_0.05_5_s0_0.01_0.mat';
files{105} = 'genLogistic_N1000_nNets100_nTrials100_maxT5k_sig0.2_sameRNG1_g1_0.05_5_s0_0.01_0.mat';

filesToPlot = [101:105];


%% For a given network, what frac of trials is it multistable

testS = 0; % s-value at which to analyze across g

tic
for ithFile = filesToPlot
    
    load([dataPath, filesep, files{ithFile}], 'numUnitsOn', 'numTrials', 'stopCriterion', 'finalStatesID', 'FIsig', 'ss', 'gs', 'N', 'numNets')
    
    % Get the index of the desired s-value
    testSind = find( (abs(ss-testS)<0.0001) );
    testNets = 1:numNets;
    
    % Calculate number of point attractor states found for each network
    netNAttSts = zeros(numel(testNets), numel(testSind), numel(gs));    %Number of point attractors at each (net, g)
    fracTrialsConvg = zeros(numNets, numel(ss), numel(gs));     % Fraction of trials that converge
    fractTrialsONatt = zeros(numNets, numel(ss), numel(gs));	% Fraction of trials converging to an active point attractor
    fractTrialsOFFatt = zeros(numNets, numel(ss), numel(gs));	% Fraction of trials converging to the quiescent point attractor
    fractTrialsPtAt2nd = zeros(numNets, numel(ss), numel(gs));	% Frac trials finding non-domainant pt att, including quiescent
    for ithNet = testNets
        for iths = 1:numel(testSind)
            for ithg = 1:numel(gs)
                % Pull data from simulation results
                netStateIDs = finalStatesID(:,ithNet,iths,ithg);
                netCondIDs = stopCriterion(:,ithNet,iths,ithg);
                isValidConvg = ismember(netCondIDs, [1,2]);
                nPointAttrs = numel(unique(netStateIDs(isValidConvg)));
                % Store the data
                netNAttSts(ithNet, iths, ithg) = nPointAttrs;
                fracTrialsConvg(ithNet, iths, ithg) = mean( isValidConvg );
                fractTrialsONatt(ithNet, iths, ithg) = mean( [isValidConvg] & [numUnitsOn(:,ithNet,iths,ithg)>0] ) ;
                fractTrialsOFFatt(ithNet, iths, ithg) = mean( [isValidConvg] & [numUnitsOn(:,ithNet,iths,ithg)==0] ) ;
                fractTrialsPtAt2nd(ithNet, iths, ithg) = mean( [isValidConvg] & [netStateIDs~=mode(netStateIDs(isValidConvg))] ) ;
            end
        end
    end
    
    X_completed = squeeze(~any(stopCriterion(:,:,testSind,:)==0, [1])); % Indicates successfully simulated g
    AttSt_count = squeeze(netNAttSts);          % Number of attractor states for each (net, g)
    multiStab_bool = squeeze(netNAttSts>1);     % Boolean indicating multistability for each (net, g)
    
    % Sort nets be center of mass of multistable g region
    tmp = multiStab_bool.*gs; tmp(tmp==0)=nan;
    [~, sortInd] = sort(nansum(tmp, 2)./sum(~isnan(tmp), 2), 'ascend');
    
    nMultiStab=sum(sum(multiStab_bool, 2)>=1);  % Number of multistable networks at each g
    
    % Plot multistability across g, Supplemental Figure 2A
    myPlotSettings(width=1.5, height=1.5)
    figure; imagesc(gs, testNets, AttSt_count(sortInd,:), 'AlphaData', X_completed(sortInd,:));
    caxis([0, 2]); colormap(parula(3)); cb.YTick=[0.33,1,1.67]; cb.YTickLabel={'0', '1', '2+'};
    xlabel('g'); ylabel('Net (sorted)')
    yticks([1, round(numNets/2), numNets])
    xticks([gs(1):1:gs(end)])
    axis square
    
    % Set up data for determining classes
    classVals = [1, 1.75, 2.6, 3.2, 3.8, 4.499, 5.25, 6, 6.5, 7]; % Color values for each class in imagesc()
    anyFailedConv = squeeze(fracTrialsConvg<1);     % All trials failed to converge
    allFailedConv = squeeze(fracTrialsConvg==0);    % Any trial failed to converge
    anyOFF = squeeze(fractTrialsOFFatt>0);          % Any trial converged to quiescent attractor
    allOFF = squeeze(fractTrialsOFFatt==1);         % All trials converged to quiescent attractor
    anyON = squeeze(fractTrialsONatt>0);            % Any trial converged to active point attractor
    allON = squeeze(fractTrialsONatt==1);           % All trials converged to active point attractor
    isMultiStab = squeeze(AttSt_count>=2);          % There were at least 2 point attractors
    isBiStab = squeeze(AttSt_count==2);             % There were only two point attractors
    
    % Determine the class for each g-value of each network
    cmapData = nan(size(allON));
    cmapData( allFailedConv                  ) = classVals(1);      % Chaos
    cmapData( anyFailedConv & anyOFF &~anyON ) = classVals(2);      % Chaos and quiescent attractor
    cmapData( allOFF                         ) = classVals(3);      % Only quiescent attractor
    cmapData( allON & ~isMultiStab           ) = classVals(4);      % Single active attractor only
    cmapData( anyFailedConv & anyON & anyOFF ) = classVals(5);      % Chaos + active attractor + quiescent + multistab
    cmapData( ~anyFailedConv & isBiStab & anyOFF&anyON &ss<(4*FIsig) ) = classVals(6);  % Single active and single quiescent attractors
    cmapData( allON &  isMultiStab           ) = classVals(7);      % Only active attractors and is multistable
    cmapData(~anyFailedConv & anyOFF&anyON&isMultiStab&~isBiStab ) = classVals(8);      % Multiple active attractors + quiescent
    cmapData( anyFailedConv & anyON & ~anyOFF & ~isMultiStab     ) = classVals(9);      % Chaos + single active attractor
    cmapData( anyFailedConv & anyON & ~anyOFF &  isMultiStab     ) = classVals(10);     % Chaos + multiple active attractors
    
    [class_counts, ~] = hist(cmapData(:),unique(classVals));
    
    % Throw a warning if any (net,g) point is an undefined class
    if any(isnan(cmapData), 'all')
        warning([num2str(sum(isnan(cmapData), 'all')), ' points'' activity regimes not defined'])
        candidateNewClass_count = sum( anyOFF&anyON&isMultiStab&anyFailedConv , 'all');
        disp(candidateNewClass_count)
    end
    
    % Plot classes across g, Supplemental Figure 2B
    myPlotSettings(width=1.5, height=1.5)
    figure; imagesc(gs, testNets, cmapData(sortInd,:), 'AlphaData', [X_completed(sortInd,:)&~isnan(cmapData)]);
    colormap(turbo); caxis([1 7])
    xlabel('g'); ylabel('Net (sorted)')
    yticks([1, round(numNets/2), numNets])
    xticks([gs(1):1:gs(end)])
    axis square
    
    % Print results for this N
    disp(['At N=', num2str(N), ':'])
    disp(['Class counts: ', num2str(class_counts)])
    disp(['Peak frac nets multistable at same g ', num2str(max(mean(multiStab_bool, 1)))])
    disp(['Frac nets multistable at any g ', num2str(nMultiStab/numNets)])
    
end
toc

% Plot a colorbar for number of point attractors
figure
cb = colorbar; clabel = 'Number of Point Att.'; ylabel(cb, clabel);
caxis([0, 2]); colormap(parula(3)); cb.YTick=[0.33,1,1.67]; cb.YTickLabel={'0', '1', '2+'};

% Plot a colorbar for classes
figure;
n=256; cm=turbo(n);
classVals_used = classVals;
cvars_map = ceil((classVals_used - min(classVals_used))/(max(classVals_used)-min(classVals_used)) * (n-1) + 1);
colormap(cm(cvars_map, :) )
cb = colorbar;
steps = (max(classVals)-min(classVals))/numel(classVals)/2;
cb.YTick= [1+1*steps:2*steps:30*steps] ; cb.YTickLabel={num2str([1:numel(classVals)]')};
caxis([min(classVals), max(classVals)])
ylabel(cb, 'Classes');

myPlotSettings
