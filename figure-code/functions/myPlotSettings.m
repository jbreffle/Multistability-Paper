function [] = myPlotSettings(varargin)
% myPlotSettings - This function sets default figure settings when creating
% a new figure. The function accepts a variable number of input
% arguments, which can be used to override the default settings.
%
%   INPUTS:
%       varargin - Optional input arguments that can be used to override the
%       default settings. The function accepts the following parameters:
%
%           width: Figure width (default: 4.25)
%           height: Figure height (default: 3)
%           lw: Line width (default: 1.5)
%           afs: Font size (default: 10)
%           ttlfsz: Title font size (default: 1.2)
%           txtfsz: Text font size (default: 10)
%           alw: Axes line width (default: 0.75)
%           msz: Marker size (default: 6)
%           afw: Axes font weight (default: normal)
%           afn: Axes font name (default: times)
%
%   OUTPUTS:
%       None - The function modifies the default plotting settings.
%
%   EXAMPLES:
%       Example 1: Override the default figure width and height
%           myPlotSettings(width=6, height=4);
%
%       Example 2: Override the default line width and font size
%           myPlotSettings('lw', 2, 'afs', 12);


%% Parse varargin
inputObj = inputParser;
inputObj.addParameter('width', 4.25);
inputObj.addParameter('height', 3);
inputObj.addParameter('lw', 1.5);
inputObj.addParameter('afs', 12);
inputObj.addParameter('ttlfsz', 1.2);
inputObj.addParameter('txtfsz', 10);
inputObj.addParameter('alw', 0.75);
inputObj.addParameter('msz', 8);
inputObj.addParameter('afw', 'Normal');
inputObj.addParameter('afn', 'times');

inputObj.parse(varargin{:});

p = inputObj.Results;

%% Set figure properties
set(groot, ...
    'DefaultLineLineWidth', p.lw, ...
    'DefaultImplicitFunctionLineLineWidth', p.lw, ...
    'DefaultErrorBarLineWidth', p.lw, ...
    'DefaultLineMarkerSize', p.msz, ...
    'DefaultAxesLineWidth', p.alw, ...
    'DefaultAxesTitleFontSize', p.ttlfsz, ...
    'defaulttextfontsize', p.txtfsz, ...
    'DefaultAxesFontSize', p.afs, ...
    'DefaultAxesFontWeight', p.afw, ... % normal, bold, 
    'defaultfigurecolor', [1 1 1], ...
    'defaultaxesFontName', p.afn, ... % times, cambria math
    'defaultFigureUnits', 'inches', ...
    'defaultFigurePaperUnits','inches', ...
    'Units', 'inches', ...
    'defaultFigurePosition', [2, 3, p.width, p.height], ...
    'defaultFigurePaperPositionMode', 'auto', ...
    'defaultFigureInvertHardCopy', 'off' ...
);
end
