% perturb_nets_run_sim
%
% Code to produce data for panel panel E of Figure 2 of the paper
% "Multistability in neural systems with random cross-connections"
%
% Limit cycle pertubation analysis for each network:
% -At 100 points in time randomly pick 10% of units and perturb their rates
%   on the order of +/- 0.0001 Hz
% -Simulate from that perturbation
% -Compare the perturbation against the original simulation:
%   sqrt(mean( (original-pertubation)^2 ) )
%
% Deviation should exponentially decrease for point attractors and
% exponentially increase for chaos.

addpath('./functions')

%% Select the directory in which to save the simulation results

saveResults = true; % If true, saves results at end of script
resultsDir = './results/'; % location to save data


%% Network parameters
N=100; s=0.0; g=3.25;    % Network parameters
isTanh = 0; FIsig = 0.1; % FI curve parameters

%% Script options
nNets = 100;            % Number of networks to simulate
nTrials = 100;          % Number of trials to simulate for each network

% Generate vector of rng seeds for each network and trial
rng('default')
netSeedVec = [1,  randi((2^32-1), 1, (nNets-1))];   % Include the example net, rng(1)
initSeedVec = [0, randi((2^32-1), 1, (nTrials-1))]; % include the 0 seed, which starts all units near zero

%% Sim parameters

% ode45 options
maxstep = [];       % Default is one-tenth of the tspan interval
RelTol = 1e-6;      % Default is 1e-3
AbsTol = RelTol/10; % Default is 1e-6

maxTime = 21 * 1000; % Duration of unperturbed sim, units of tau

maxUnitChange_tau = 10; % In the final maxUnitChange_tau of the sim, what is the max range of any unit's rate?

%% Perturbation parameters
fracPerturb = 0.1;          % Fraction of units to perturb
perturbMagnitude = 0.00001; % Amount perturbed units are perturbed +/- in Hz
nPerturbations = 100;       % Number of perturbations per trial
pertMaxT = 200; partTStepSize = 0.001; % Perturbation time vec params for averaging over
minPertT = (20/21)*maxTime;            % Don't perturb before this time point
perturbTimes = linspace(minPertT, maxTime-pertMaxT, nPerturbations+1); perturbTimes = perturbTimes(1:end-1);
nPerturbed = max( round(N*fracPerturb) , 1); % Number of perturbed units

linfitend = 200; % Range of time (tau) following perturbation to fitlm to
nRMSEBins = linfitend*10;
RMSrateEndAvg = linfitend; % taus from end over which to average for RMSAvgFinalRates

%% Automatic set up

% Initial RMSE value
expectedRMSE = sqrt(mean ( [repmat(perturbMagnitude, 1, nPerturbed), zeros(1, N-nPerturbed)].^2));

% Set up logistic functions
logisticThresh  = 0.5 + sqrt(0.25-FIsig) + FIsig*log(FIsig)-2*FIsig*log(0.5 + sqrt(0.25-FIsig));
logisticFun = @(y) (1./(1+exp( -(y-logisticThresh)/FIsig))); % I=input -> f(I)=rate
invLogisticFun = @(r) logisticThresh - FIsig*log( (1./r) - 1) ; % F(I)=rate -> I=input

% Set up tanh functions
tanhThresh =  -1* (sqrt(1-FIsig) -0.5*FIsig*log( (1 + sqrt(1-FIsig))/(1-sqrt(1-FIsig)) ));
tanhFun = @(y) (tanh( (y-tanhThresh) / FIsig ) );
invTanhFun = @(r) tanhThresh + FIsig*atanh(r);

% Parameters for random initialization of units' rates
init_mu = 0.5; init_alpha = 0.1;
initSigmoid = @(y) (1./(1+exp( -(y-init_mu)/init_alpha)));


%% Main loop

% Initialze data matrices
fitR2Mat_fixedb = zeros(nNets, nTrials);        % R2 of fit to RMSdeviation with fixed intercept
fitSlopetMat_fixedb = zeros(nNets, nTrials);    % Slope of fit to RMSdeviation with fixed intercept
fitR2Mat = zeros(nNets, nTrials);               % R2 of fit to RMSdeviation without fixed intercept
fitSlopetMat = zeros(nNets, nTrials);           % Slope of fit to RMSdeviation without fixed intercept
meanRMSEMat = zeros(nNets, nTrials, nRMSEBins); % matrix of the median RMSdeviation for each trial of each ent
maxUnitVal_y = zeros(nNets, nTrials);           % Maximum rate of any unit at the end of the simulation
maxUnitChange_y = zeros(nNets, nTrials);        % Maximum change in any units firing rate in the last maxUnitChange_tau of the simulation
RMSFinalRates = zeros(nNets, nTrials);          % RMS of final rate
RMSAvgFinalRates = zeros(nNets, nTrials);       % RMS of average of final rate
maxUnitVal_perts = zeros(nNets, nTrials, nPerturbations); % Maximum rate of any unit at the end of each perturbation

% Loop over networks
globalStartTime=tic;
for ithNet = 1:numel(netSeedVec)
    disp(['Starting net ', num2str(ithNet), ', current total runTime=', num2str(toc(globalStartTime))])
    netSeed = netSeedVec(ithNet);
    
    %% Create ithNet network
    rng(netSeed, 'twister');    % the state of the rng for initializing W
    W = sqrt(1/N)*randn(N);     % variance = sigma^2 = 1/N, std = sqrt(variance) = sqrt(1/N)
    W = W - diag(diag(W));
    W = W';
    
    % Loop over trials
    parfor ithTrial = 1:numel(initSeedVec)
    %for ithTrial = 1:numel(initSeedVec)
        
        %% Set up for chosen FI curve
        if isTanh
            transFun = tanhFun;
            invTransFun = invTanhFun;
            rateBounds = [-1+eps, 1-eps];
        else
            transFun = logisticFun;
            invTransFun = invLogisticFun;
            rateBounds = [0+eps, 1-eps];
        end

        %% Set initial condtions
        initSeed = initSeedVec(ithTrial);
        rng(initSeed, 'twister');  % the state of the rng for the initial condition
        
        r_init = (initSigmoid(rand(N,1))-0.5)*diff(rateBounds)+mean(rateBounds);
        I_init = invTransFun( r_init );
        
        if initSeed==0
            I_init = zeros(1,N);
        end
        
        %% Run unperturbed simulation
        odefun = @(t, y) ( -y + s*transFun(y) + g*(W*transFun(y)) );
        odeOptions = odeset('MaxStep', maxstep, 'RelTol',RelTol,'AbsTol',AbsTol);
        [t,y] = ode45(odefun, [0, maxTime], I_init, odeOptions);
        
        % Store results
        unitsMaxChange = max(transFun(y([t>(maxTime-maxUnitChange_tau)],:)))-min(transFun(y([t>(maxTime-maxUnitChange_tau)],:)));
        maxUnitChange_y(ithNet, ithTrial) = max(abs(unitsMaxChange));
        maxUnitVal_y(ithNet, ithTrial) = transFun(max(y(end,:)));
        RMSFinalRates(ithNet, ithTrial) = sqrt( mean( transFun(y(end,:)).^2 ));
        RMSAvgFinalRates(ithNet, ithTrial) = sqrt( mean( mean(transFun(y([t>(maxTime-RMSrateEndAvg)],:))).^2, 2)) ;
        
        %% Perturbation analysis
        interp_t = [0:partTStepSize:pertMaxT]; % Time vector of points at which to interpolate ode45 results
        deviationOP = nan( numel(perturbTimes), numel(interp_t) );  % Matrix of RMSdeviations across time for each perturbation
        trial_maxUnitVal_perts = zeros(1, nPerturbations);
        for ithPert = 1:numel(perturbTimes)
            
            % Select units and the time index for ithPert perturb
            perturbInds = datasample( [1:N], nPerturbed );
            [~, startInd] = min(abs(t-  perturbTimes(ithPert)           ));
            [~, endInd] =   min(abs(t- (perturbTimes(ithPert)+pertMaxT) )); endInd = endInd+1;
            t_p = t(startInd:endInd);
            
            % Perturb initial condition rates
            R_init_p = transFun(y(startInd,:));
            R_perturb_vec = perturbMagnitude .* (2*[rand(1,nPerturbed)<0.5]-1);
            % If a perturbation takes a unit's rate outside its rate bounds, then flip the sign of that perturbation
            outOfBoundsInds = (R_init_p(perturbInds)+R_perturb_vec)<rateBounds(1) | (R_init_p(perturbInds)+R_perturb_vec)>rateBounds(2);
            R_perturb_vec(outOfBoundsInds) = -1*R_perturb_vec(outOfBoundsInds);
            R_init_p(perturbInds) = R_init_p(perturbInds) + R_perturb_vec;
            R_init_p(R_init_p>rateBounds(2)) = rateBounds(2)-eps; R_init_p(R_init_p<rateBounds(1))=rateBounds(1)+eps;
            I_init_p = invTransFun(R_init_p);
            
            % Simulate the perturbation
            [~,y_p] = ode45(odefun, t_p, I_init_p, odeOptions);
            trial_maxUnitVal_perts(ithPert) = transFun(max(y_p(end,:)));
            
            % Append to the simulation results, to get correct matrix size
            if size(y,1)-size(y_p,1)-startInd+1>0
                y_p = [y_p; nan(size(y,1)-size(y_p,1)-startInd+1, N)];
            end
            op = [y(1:startInd-1,:); y_p];
            
            % Calculate interpolated deviation from original simulation rates
            y_diff = transFun(y_p) - transFun(y(startInd:end,:));
            dev = sqrt(mean( (y_diff).^2, 2));
            deviationOP(ithPert,:) = interp1(t_p-t_p(1), dev(1:numel(t_p)), interp_t );
            
        end % Perturbation loop
        
        %% Calculate RMSdeviation linear fit
        startInd = 1; [~, endInd] = min(abs(interp_t-linfitend));
        x = interp_t(startInd:endInd);
        yToFit = log(median(deviationOP(:,startInd:endInd), 1));
        include=~isnan(yToFit) & ~isinf(yToFit);
        
        % Calculate correlation with fixed b value (the intercept)
        [mdl,gof,~] = fit(x(include)', yToFit(include)', 'poly1', 'Lower',[-inf, log(expectedRMSE) ],'Upper',[inf, log(expectedRMSE) ]);
        fitCoeffs_fixedb = coeffvalues(mdl);
        mdlRsqr_fixedb = gof.rsquare;
        % Calculate correlation without fixed b value (the intercept)
        mdl = fitlm(x(include), yToFit(include));
        fitCoeffs = fliplr(mdl.Coefficients.Estimate');
        mdlRsqr = mdl.Rsquared.Ordinary;
        
        edges = x(1:round(endInd/nRMSEBins):endInd);
        [~,~,loc]=histcounts(x,edges);
        meany = accumarray(loc(:),yToFit(:))./accumarray(loc(:),1);
        xmid = 0.5*(edges(1:end-1)+edges(2:end));
        
        % Store results
        maxUnitVal_perts(ithNet, ithTrial, :) = trial_maxUnitVal_perts;
        fitR2Mat_fixedb(ithNet, ithTrial) = mdlRsqr_fixedb;
        fitSlopetMat_fixedb(ithNet, ithTrial) = fitCoeffs_fixedb(1);
        fitR2Mat(ithNet, ithTrial) = mdlRsqr;
        fitSlopetMat(ithNet, ithTrial) = fitCoeffs(1);
        meanRMSEMat(ithNet, ithTrial, :) = exp(meany');
        
    end % initSeed loop
end % netSeed loop
globalRunTime = toc(globalStartTime);
disp(['Global run time: ', num2str(datestr(datenum(0,0,0,0,0,globalRunTime),'HH:MM:SS')), ' (HH:MM:SS)'])

%% Save results

% Some calculations of the median RMSE vals for convenience
initRMSE = mean(meanRMSEMat(:,:,1:2), 3);       % Mean of first two bins
endRMSE = mean(meanRMSEMat(:,:,end-1:end), 3);  % Mean of last two bins
diffRMSE = endRMSE-initRMSE;
maxRMSE = max(meanRMSEMat, [], [3]);

if saveResults
    saveTime = char(datetime('now', 'format', 'yyyy-MM-dd'));
    saveName = [saveTime, ...
        '_nNets', num2str(nNets), ...
        '_nTrials', num2str(nTrials)...
        '_maxT', num2str(round(maxTime/1000)), 'k', ...
        '_log10pertMag', num2str(round(log10(perturbMagnitude))), ...
        '_s', num2str(s), ...
        '_g', num2str(g), ...
        '_FIsig', num2str(FIsig), ...
        '_pert.mat'];
    savePath = [resultsDir filesep saveName];
    safeSavePath = AppendFileNum(savePath);
    
    clear op y deviationOP mdl y_diff y_p interp_t loc x yToFit t include
    save(safeSavePath)
else
    warning('Did not save results')
end
