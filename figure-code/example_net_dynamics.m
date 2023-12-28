% example_net_dynamics.m
%
% Code to produce either 
%   1) panels A-D of Figure 2 or 
%   2) panel A of Figure 3 of the paper
% "Multistability in neural systems with random cross-connections"
%
% Limit cycle pertubation analysis:
% -At 100 points in time randomly pick 10% of units and perturb their rates
%   on the order of +/- 0.0001 Hz
% -Simulate from that perturbation
% -Compare the perturbation against the original simulation:
%   sqrt(mean( (original-pertubation)^2 ) )
%
% Deviation should exponentially decrease for point attractors and
% exponentially increase for chaos.

addpath('./functions')

%% Select Figure to plot

figSim = 'fig3'; % Select 'fig2' or 'fig3'


%% Simulation parameters

switch figSim
    case 'fig2' % Plots Figure 2A-D
        % Network with all types of activity
        N=100; s=0.0; g=3.25;    % Network parameters
        isTanh = 0; FIsig = 0.1; % FI curve parameters
        netSeed = 1;             % Seed for network RNG
        % initSeedVec: 0 for quiescent; 1,4 for PtAtt; 2 for chaos; 7,9 for limit cycle
        initSeedVec = [0,1,4,2,7,9]; % RNG seeds used for initial conditions
        runPerturb      = 1; % Flag to choose to run perturbations
        runAppendedPCA  = 1; % Flag to choose to run trial-appended PCA
    case 'fig3' % Plots Figure 3A
        % Example of net with only on and off attractor state (checked first 100 trials)
        N=100; s=0.0; g=2.0;     % Network parameters
        isTanh = 0; FIsig = 0.1; % FI curve parameters
        netSeed = 5;             % Seed for network RNG
        initSeedVec = [0:1]; % RNG seeds used for initial conditions
        runPerturb      = 0; % Flag to choose to run perturbations
        runAppendedPCA  = 0; % Flag to choose to run trial-appended PCA
    otherwise
        error('Unkonwn case')
end

maxTime = 21 * 1000; % Simulation duration for unperturbed simulation

% ode45 options
maxstep = 1;
RelTol = []; % Default is 1e-3
AbsTol = []; % Default is 1e-6

% Perturbation parameters
fracPerturb = 0.1;          % Fraction of units to perturb
perturbMagnitude = 0.00001; % Amount perturbed units are perturbed +/- in Hz
nPerturbations = 100;       % Number of perturbations per trial
pertMaxT = 200;             % Duration of perturbation to simulate
partTStepSize = 0.001;      % Timestep for interpolating perturbation magnitude
minPertT = (20/21)*maxTime; % Don't perturb before this time point
perturbTimes = linspace(minPertT, maxTime-pertMaxT, nPerturbations+1); perturbTimes = perturbTimes(1:end-1); % Vector of perturbation start times
nPerturbed = max( round(N*fracPerturb) , 1); % Number of perturbed units


%% Automatic set up

% Initial RMSE value
expectedRMSE = sqrt(mean ( [repmat(perturbMagnitude, 1, nPerturbed), zeros(1, N-nPerturbed)].^2));

% Set up logistic functions
logisticThresh  = 0.5 + sqrt(0.25-FIsig) + FIsig*log(FIsig)-2*FIsig*log(0.5 + sqrt(0.25-FIsig));
logisticFun = @(y) (1./(1+exp( -(y-logisticThresh)/FIsig)));   % I=input -> f(I)=rate
invLogisticFun = @(r) logisticThresh - FIsig*log( (1./r) - 1); % F(I)=rate -> I=input

% Set up tanh functions
tanhThresh =  -1* (sqrt(1-FIsig) -0.5*FIsig*log( (1 + sqrt(1-FIsig))/(1-sqrt(1-FIsig)) ));
tanhFun = @(y) (tanh( (y-tanhThresh) / FIsig ) );
invTanhFun = @(r) tanhThresh + FIsig*atanh(r);

if isTanh
    transFun = tanhFun;
    invTransFun = invTanhFun;
    rateBounds = [-1+eps, 1-eps];
else
    transFun = logisticFun;
    invTransFun = invLogisticFun;
    rateBounds = [0+eps, 1-eps];
end

% Parameters for random initialization of units' rates
init_mu = 0.5; init_alpha = 0.1;
initSigmoid = @(y) (1./(1+exp( -(y-init_mu)/init_alpha)));


%% Generate the newtork
rng(netSeed);  % The state of the rng for initializing W
W = sqrt(1/N)*randn(N);
W = W - diag(diag(W));
W = W';                                    % Connection matrix without s and g
W_all = W*g + diag(diag(ones(size(W))))*s; % Connection matrix including s and g

if isequal(figSim, 'fig2')
    % Calculate example net statistics
    ccidx = ~logical(eye(size(W)));
    ccs = W(ccidx); % All cross-connections only
    x = min(ccs): (max(ccs)-min(ccs))/numel(ccs) :max(ccs);
    normDist = makedist('Normal' ,'mu',0,'sigma', 1/sqrt(N));
    normCDF = cdf( normDist,x);
    [H,P,KSSTAT,CV] = kstest(ccs, 'CDF', [x', normCDF']);
    netEigs = eig(W_all);
    disp(['Cross-connection KS-test: pval=', num2str(P), ' KS-stat=', num2str(KSSTAT)])

    myPlotSettings(width=2.25, height=1.5)

    % cc distribution, Supplemental Figure 1A
    figure; histogram(ccs);
    box off
    xlabel('Cross-connection strength'); ylabel('Cross-connections (count)')
    disp(['CC Mean: ', num2str(mean(ccs)), ', STD: ', num2str(std(ccs))])

    % Eigenvalue scatter plot, Supplemental Figure 1B
    figure; hold on;
    viscircles([s, 0], 1*g, 'Color', 'k');
    scatter(real(netEigs), imag(netEigs))
    axis equal
    xlabel('Real part'); ylabel('Imaginary part')
    myPlotSettings(width=3, height=2)
end

%% Main Loop

% Data for appended PCA. All simulated trials are plotted in the same PCA space
appendedPCA_y = [];
appendedPCA_t = [];
appendedPCA_i = [];

globalTic = tic;
for initSeed = initSeedVec
    
    %% Set initial condtions
    rng(max(initSeed, 0));
    
    r_init = (initSigmoid(rand(N,1))-0.5)*diff(rateBounds)+mean(rateBounds);
    I_init = invTransFun( r_init );
    
    if initSeed==-1
        I_init = zeros(1,N);
    end
    if initSeed==0
        I_init = 1/2*randn(1,N);
    end
    
    %% Unperturbed simulation
    odefun = @(t, y) ( -y + s*transFun(y) + g*(W*transFun(y)) );
    odeOptions = odeset('MaxStep', maxstep, 'RelTol',RelTol, 'AbsTol',AbsTol);
    
    tic
    [t,y] = ode45(odefun, [0, maxTime], I_init, odeOptions);
    simTime = toc;
    disp(['Simulation duration was ', num2str(simTime), ' seconds'])
    
    if runAppendedPCA
        appendedPCA_y = [appendedPCA_y; y];
        appendedPCA_t = [appendedPCA_t; t];
        appendedPCA_i = [appendedPCA_i; repmat(find(initSeed==initSeedVec), size(t))];
    end
    
    %% Plot unperturbed simulation
    myPlotSettings(width=1.5, height=1)
    
    rateOffset = 0.05; % yaxis expansion of rateBounds
    unitsToPlot = 1:min(25, N); % Subset of the N units to plot
    
    % Plot sections of time, with random subset of units
    figure(1000+initSeed);
    
    % Choose how to plot based on known result of each trial
    if isequal(figSim, 'fig3')
        timeWindow = 30;
        startInd = 1;
        [~,endInd] = min(abs(t-(0+timeWindow)));
        timeToPlot = [startInd:endInd];
        plot(t(timeToPlot), transFun(y(timeToPlot, unitsToPlot)))
        ylim([rateBounds+[-rateOffset, rateOffset]]); xlim([0, timeWindow])
        xticks([0, 15, 30])
    elseif ismember(initSeed, [0,1,4]) % Figure 2 point attractors
        timeWindow = 30;
        startInd = 1;
        [~,endInd] = min(abs(t-(0+timeWindow)));
        timeToPlot = [startInd:endInd];
        plot(t(timeToPlot), transFun(y(timeToPlot, unitsToPlot)))
        ylim([rateBounds+[-rateOffset, rateOffset]]); xlim([0, timeWindow])
        xticks([0, 15, 30])
    else % Figure 2 limit cycle or chaos
        timeWindow = 100;
        [~,startInd] = min(abs(t-(maxTime-timeWindow)));
        [~,endInd] = min(abs(t-maxTime));
        timeToPlot = [startInd:endInd];
        plot(t(timeToPlot), transFun(y(timeToPlot, unitsToPlot)))
        ylim([rateBounds+[-rateOffset, rateOffset]]); xlim([maxTime-timeWindow, maxTime])
        xticks([maxTime-timeWindow, maxTime])
    end
    xlabel('time (\tau)'); ylabel('f');
    yticks([rateBounds])
    ax = gca;
    ax.XAxis.TickLength = [0 0];
    box off
    myPlotSettings(width=3, height=2)
    
    %% Perturbation analysis
    if runPerturb
        tic
        tInterp = [0:partTStepSize:pertMaxT]; % Time vector of points at which to interpolate ode45 results
        deviationOP = nan( numel(perturbTimes), numel(tInterp) ); % Matrix of RMSdeviations across time for each perturbation
        for ithPert = 1:numel(perturbTimes)
            
            % Select units and the time index for ithPert perturb
            perturbInds = datasample( [1:N], nPerturbed );
            [~, startInd] = min(abs(t-  perturbTimes(ithPert)           ));
            [~, endInd] =   min(abs(t- (perturbTimes(ithPert)+pertMaxT) )); endInd = endInd+1;
            tPerturbation = t(startInd:endInd);
            
            % Perturb initial condition rates
            RInitPerturbation = transFun(y(startInd,:));
            R_perturb_vec = perturbMagnitude .* (2*[rand(1,nPerturbed)<0.5]-1);
            % If a perturbation takes a unit's rate outside its rate bounds, then flip the sign of that perturbation
            outOfBoundsInds = (RInitPerturbation(perturbInds)+R_perturb_vec)<rateBounds(1) | (RInitPerturbation(perturbInds)+R_perturb_vec)>rateBounds(2);
            R_perturb_vec(outOfBoundsInds) = -1*R_perturb_vec(outOfBoundsInds);
            RInitPerturbation(perturbInds) = RInitPerturbation(perturbInds) + R_perturb_vec;
            RInitPerturbation(RInitPerturbation>rateBounds(2)) = rateBounds(2)-eps; RInitPerturbation(RInitPerturbation<rateBounds(1))=rateBounds(1)+eps;
            I_init_p = invTransFun(RInitPerturbation);
            
            % Simulate the perturbation
            [~,yPerturbation] = ode45(odefun, tPerturbation, I_init_p, odeOptions);
            
            % Append to the simulation results, to get correct matrix size
            if size(y,1)-size(yPerturbation,1)-startInd+1>0
                yPerturbation = [yPerturbation; nan(size(y,1)-size(yPerturbation,1)-startInd+1, N)];
            end
            
            % Calculate interpolated deviation from original simulation rates
            y_diff = transFun(yPerturbation) - transFun(y(startInd:end,:));
            dev = sqrt(mean( (y_diff).^2, 2));
            deviationOP(ithPert,:) = interp1(tPerturbation-tPerturbation(1), dev(1:numel(tPerturbation)), tInterp );
        end
        
        %% Calculate RMSdeviation linear fit
        linfitend = 100; % Range of time (tau) following perturbation to fitlm to
        startInd = 1; [~, endInd] = min(abs(tInterp-linfitend));
        x = tInterp(startInd:endInd);
        yToFit = log10(median(deviationOP(:,startInd:endInd), 1));
        
        include=~isnan(yToFit) & ~isinf(yToFit);
        [mdl,gof,~] = fit(x(include)', yToFit(include)', 'poly1', 'Lower',[-inf, log10(expectedRMSE) ],'Upper',[inf, log10(expectedRMSE) ]);
        fitCoeffs = coeffvalues(mdl);
        mdlRsqr = gof.rsquare;
        
        %% RMSE plot, Figure 2D
        myPlotSettings(width=2.95, height=1.5)
        cm = colormap(turbo(numel(unique(initSeedVec))));
        figure(9000); hold on;
        tmp_plot = plot(tInterp, (median(deviationOP, 1)), 'Color', cm(find(initSeedVec==initSeed),:));
        if initSeed~=2 % Keep the chaotic trial plot on top
            uistack(tmp_plot,'bottom')
        end
        yline( (expectedRMSE), ':' )
        xlabel('Time from perturbation (\tau)'); ylabel('RMS deviation')
        set(gca,'YScale', 'Log')
        xlim([0, 100])
        ylim([(eps/50), 1])
        myPlotSettings(width=3, height=2)
        
        pertRunTime = toc;
        disp(['Perturbation run time: ', num2str(pertRunTime)])
    end
    
end % initSeed loop
globalRunTime = toc(globalTic)


%% Run PCA on all completed trials
if runAppendedPCA
    myPlotSettings(width=2.75, height=1.3)
    
    % Choose time range to plot
    appendedPCA_timeSubset = true;
    PCA_minT = maxTime*0.95;
    PCA_maxT = maxTime*1.0;
    if appendedPCA_timeSubset
        validInds = appendedPCA_t>=PCA_minT & appendedPCA_t<=PCA_maxT;
    else
        validInds = 1:numel(appendedPCA_t);
    end
    
    % Calculate PCA for choosen (sub)set of data
    dataToPlot = transFun(appendedPCA_y(validInds,:));
    colorVec = appendedPCA_i(validInds)';
    [~,PCAdataToPlot, ~, TSQUARED, EXPLAINED, ~] = pca(dataToPlot);
    
    % Plot PCA of all trials
    figure; hold on;
    view([-75, 25]);
    cm = colormap(turbo(numel(unique(initSeedVec))));
    xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
    for ii = 1:numel(initSeedVec)
        trialInds = [appendedPCA_i(validInds)==ii];
        
        trialStart = find(trialInds, 1, 'first');
        trialEndInd = find(trialInds, 1, 'last');
        plot3(PCAdataToPlot(trialInds,1), PCAdataToPlot(trialInds,2), PCAdataToPlot(trialInds,3), 'Color', cm(ii,:))
        
        trialInds = find(trialInds);
        if diff(PCAdataToPlot(trialInds(end-100:end),1))<0.001
            scatter3(PCAdataToPlot(trialEndInd,1), PCAdataToPlot(trialEndInd,2), PCAdataToPlot(trialEndInd,3), 35,  cm(ii,:), 'o', 'filled')
        end
    end
    grid on
    
    % Plot the chaotic trial, Supplemental Figure 1C
    myPlotSettings(width=2.25, height=1.5)
    chaoticTrialSeed = 2;
    chaoticTrialInds = [appendedPCA_i==find( initSeedVec==chaoticTrialSeed )];
    [~,chaoticPCA, ~, ~, ~, ~] = pca(transFun(appendedPCA_y(chaoticTrialInds,:)));
    
    figure; hold on; view(-20, 20)
    xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
    patch([chaoticPCA(:,1)' nan], [chaoticPCA(:,2)' nan], [chaoticPCA(:,3)' nan], [appendedPCA_t(chaoticTrialInds)' nan], 'EdgeColor','interp','FaceColor','none')
    colormap(turbo)
    grid on
    
    % Stand-alone colorbor corresponding to Supplemental Figure 1C
    figure;
    colormap(turbo)
    cb = colorbar(); ylabel(cb, 'time (\tau)');
    caxis([0, maxTime])
    cb.Ticks = [0, maxTime];
    cb.TickLabels={cb.Ticks(1), ['10^', num2str(log10(cb.Ticks(end)))]};
    
    % Calculate mean and STD of RMSrates
    minT = maxTime*0.95;
    finalMeanRates = zeros(size(initSeedVec));
    finalMeanRatesSTD = zeros(size(initSeedVec));
    for ii = 1:numel(initSeedVec)
        trialInds = [appendedPCA_i==ii];
        trialRates = transFun(appendedPCA_y(trialInds&appendedPCA_t>minT,:));
        avgNetTrialRates = sqrt(mean(trialRates.^2, 2));
        finalMeanRates(ii) = mean(avgNetTrialRates);
        finalMeanRatesSTD(ii) = std(avgNetTrialRates);
    end
    
    % Plot RMSrates, Figure 2C
    myPlotSettings(width=2.75, height=1.3)
    figure; hold on
    for ii = 1:numel(initSeedVec)
        errorbar(ii, finalMeanRates(ii), finalMeanRatesSTD(ii), 'o', 'CapSize', 20, 'Color', cm(ii,:), 'MarkerFaceColor', cm(ii,:), 'LineWidth', 2)
    end
    xticks(1:6); xticklabels({'OFF', 'ON', 'ON', 'Chaos', 'Cycle', 'Cycle'})
    xticks(1:6); xticklabels({'Quiescent', 'Stable active', 'Stable active', 'Chaos', 'Cycle', 'Cycle'})
    xlim([0.6, 6.4])
    ylim([-0.0, 0.65]); box off; yticks([0, 0.5])
    ylabel('RMS rates')
    myPlotSettings(width=3, height=2)
    
end
