% test_zipf_plot.m
%
% Code plots data on number of visits to each fixed point on a log-log
% scale to assess whether there is a power law.
% x-axis is ranking of attractors in descending number of visits
% y-axis is number of visits of that attractor
% If an attractor visit were equivalent to a word production then this
% would be a manifestation of Zipf's Law.
%
% Paul Miller, June 5, 2023. pmiller@brandeis.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
load('bin4a_net3_N200_connflag_1_s0.5_g1.5_rnd301_starts10000000.mat')

sort_visits = sort(visits','descend');
ls = log10(sort_visits);            % approx size of attractors
ln = log10(1:length(ls))';          % rank of attractor size

% To avoid the fitted line being dominated by the few states with
% (exponentially) largest nymber of visits but fit the bulk of the data to
% a straight line better:
fitvals = find(sort_visits>5)';
num_fit = length(fitvals);

[f g] = fit(ln(fitvals),ls(fitvals),'poly1')


set(0,'DefaultLineLineWidth', 3);
set(0,'DefaultAxesFontSize',20);

set(0,'DefaultLineMarkerSize',10);



figure(3)
clf
hold on
plot(ln,f.p2+ln*f.p1,'--')      % fitted straight line
plot(ln,ls,'.')                 % actual data 

set(gca,'XTick',[0 1 2 3 4])
set(gca,'YTick',[0 2 4 6])
set(gca,'XTickLabel',{'1', '10', '100', '1000', '10^{4}'})
set(gca,'YTickLabel',{'1', '100', '10^{4}', '10^{6}'})

axis([0 4.5 -0.5 4])
xlabel('Rank')
ylabel('No. of Initial Conditions')


