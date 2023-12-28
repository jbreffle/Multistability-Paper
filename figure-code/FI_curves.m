% FI_curves.m
%
% Code to produce Figure 1 of the paper "Multistability
% in neural systems with random cross-connections"
%
% Plots logistic, tanh, and binary FI curves.

addpath('./functions')

%% Set up FI curve parameters and functions

% FI curve parameters
logisticDelta = 0.1;
tanhDelta = 1;
binaryThresh = 1;

% Calculate logisticThresh such that s=1 is the threshold for single-unit bistability
syms x s mu
f = -x +  1/(1+exp( -(x-mu) /logisticDelta ));
dfdx = diff(f, x);   % compute derivative to find critical points
IssCriticalPoints = vpasolve(f,dfdx, 0.7);
logisticThresh = IssCriticalPoints.mu;

% Calculate tanhThresh such that s=1 is the threshold for single-unit bistability
syms x s mu
f = -x + tanh( (x-mu) /tanhDelta );
dfdx = diff(f, x);   % compute derivative to find critical points
IssCriticalPoints = vpasolve(f,dfdx);
tanhThresh = IssCriticalPoints.mu;

% Define the three FI functions
logisticFun = @(y) (1./(1+exp( -(y-logisticThresh)/logisticDelta)));
tanhFun = @(y) tanh( (y-tanhThresh)/tanhDelta);
binaryFun = @(y) y>binaryThresh;

%% Plot the three panels of Figure 1
myPlotSettings(width=1.9, heigh=1.5, lw=2)

xVec = -2:0.025:2;

% Plot logistic FI
figure; hold on
plot(0.5*xVec, xVec, '--'); plot(1.0*xVec, xVec, '--'); plot(1.5*xVec, xVec, '--');
plot(xVec, logisticFun(xVec));
xlim([-0.5, 1.5]); xticks([-0.5:0.5:1.5])
ylim([-0.1, 1.1]); yticks([0, 0.5, 1])
xlabel('x'); ylabel('f(x)');
box on

% Plot tanh FI
figure; hold on
plot(0.5*xVec, xVec, '--'); plot(1.0*xVec, xVec, '--'); plot(1.5*xVec, xVec, '--');
plot(xVec, tanhFun(xVec));
xlim([-1.5, 1.5]); xticks([-1, 0, 1])
ylim([-1, 1]); yticks([-1, 0, 1])
xlabel('x'); ylabel('f(x)');
box on

% Plot binary FI
figure; hold on
plot(0.5*xVec, xVec, '--'); plot(1.0*xVec, xVec, '--'); plot(1.5*xVec, xVec, '--');
plot(xVec, binaryFun(xVec));
xlim([-0.5, 1.5]); xticks([-0.5:0.5:1.5])
ylim([-0.1, 1.1]); yticks([0, 0.5, 1])
xlabel('x'); ylabel('f(x)');
box on

% Return plot settings to default
myPlotSettings
