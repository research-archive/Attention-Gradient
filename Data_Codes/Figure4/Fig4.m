
%% ========================================================================
%  Figure 4 – EJN (Revision 1)
% -------------------------------------------------------------------------
%  This script generates all panels and statistical analyses associated
%  with Figure 4 of the manuscript.
%
%  Panels include:
%    Panel B – Correlation between cortical thickness 
%               (Left IFG + PreCG) and attentional gradient
%    Panel D – Thickness × Trait Anxiety interaction (median split)
%    Gender-based moderation analyses
%    Panel E – Prediction plot with 95% prediction interval
%    Panel F – Permutation testing histogram
%
%  Required input files (must be in current directory):
%   Fig4_1.mat
%   
%
%   Figure 4 panels (B, D, E, F); 300 dpi 
%   Statistical values printed to MATLAB Command Window
%
%  Date: 24 Feb 2026
% ========================================================================

%% Clear workspace 
clearvars;
close all;
clc;

%% ========================================================================
%% PANEL B — Scatter Plot: Left (IFG + PreCG) Thickness
% ========================================================================

load('Fig4_1.mat');

fg1 = figure('Position', [0 0 400 400]);
hold on;

scatter( ...
    zscore(thickness), ...
    zscore(CompositeGradientFig4), ...
    20*4.5, ...
    'MarkerEdgeColor', [0.5 0.5 0.5], ...
    'MarkerFaceColor', [0.2 0.2 0.2], ...
    'LineWidth', 1.25);

% Least-squares regression line
ls = lsline;
ls.LineStyle = ':';
ls.LineWidth = 4.75;
ls.XData     = [-2.3 2.3];
ls.Color     = [0 0 0];

xlabel('Cortical thickness (z score)', 'FontSize', 14);
ylabel('Att gradient (z score)', 'FontSize', 14);

set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'FontSize', 14);
xlim([-2.5 2.5]);
ylim([-2.5 2.5]);
axis square;
hold off;
exportgraphics(fg1, 'Figure4B.tif', 'Resolution', 300);

%% ---------- Correlation statistics (reported in manuscript)

X  = thickness;
Y  = CompositeGradientFig4;
n  = length(X);
df = n - 2;

[r, p] = corr(X, Y, 'Type', 'Pearson');

% 95% CI (Fisher transformation)
z      = atanh(r);
se     = 1 / sqrt(n - 3);
z_crit = 1.96;
z_CI   = [z - z_crit*se, z + z_crit*se];
r_CI   = tanh(z_CI);

fprintf('\n=== Left IFG + PreCG Thickness ===\n');
fprintf('r = %.4f\n', r);
fprintf('df = %d\n', df);
fprintf('95%% CI = [%.4f, %.4f]\n', r_CI(1), r_CI(2));
fprintf('p = %.4f\n\n', p);

clear X Y

%% ========================================================================
%% PANEL D — Thickness × Trait Anxiety Interaction (Median Split)
% ========================================================================

idxHA = TraitAnxietyFig4 >= median(TraitAnxietyFig4);
idxLA = TraitAnxietyFig4 <  median(TraitAnxietyFig4);

fg2 = figure('Position', [0 0 400 400]);
hold on;

% High Anxiety
scatter( ...
    zscore(thickness(idxHA)), ...
    zscore(CompositeGradientFig4(idxHA)), ...
    20*4.5, ...
    'MarkerEdgeColor', [0.5 0.5 0.5], ...
    'MarkerFaceColor', [1 0 0], ...
    'LineWidth', 1.25);

% Low Anxiety
scatter( ...
    zscore(thickness(idxLA)), ...
    zscore(CompositeGradientFig4(idxLA)), ...
    20*4.5, ...
    'MarkerEdgeColor', [0.5 0.5 0.5], ...
    'MarkerFaceColor', [0 0 1], ...
    'LineWidth', 1.25);

% Regression lines
ls = lsline;
ls(1).LineStyle = ':'; ls(1).LineWidth = 4.75; ls(1).XData = [-2.3 2.3]; ls(1).Color = [0 0 1];
ls(2).LineStyle = ':'; ls(2).LineWidth = 4.75; ls(2).XData = [-2.3 2.3]; ls(2).Color = [1 0 0];

xlabel('rGMV (z score)');
ylabel('Att gradient (z score)');
set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'FontSize', 14, 'XTick', -2:1:2);

xlim([-2.5 2.5]);
ylim([-2.5 2.5]);
axis square;
hold off;
exportgraphics(fg2, 'Figure4D.tif', 'Resolution', 300);
%% ---------- Slopes and correlations for interaction

X1 = thickness(idxHA); Y1 = CompositeGradientFig4(idxHA);
X2 = thickness(idxLA); Y2 = CompositeGradientFig4(idxLA);

coef1 = polyfit(X1, Y1, 1);
coef2 = polyfit(X2, Y2, 1);

fprintf('Slope (High Anxiety) = %.4f\n', coef1(1));
fprintf('Intercept (High Anxiety) = %.4f\n', coef1(2));
fprintf('Slope (Low Anxiety) = %.4f\n', coef2(1));
fprintf('Intercept (Low Anxiety) = %.4f\n\n', coef2(2));

%% ========================================================================
%% Gender-Based ANCOVA
% ========================================================================

idxMale   = genderFig4 == 1;
idxFemale = genderFig4 == 2;

X_m = thickness(idxMale);
Y_m = CompositeGradientFig4(idxMale);

X_f = thickness(idxFemale);
Y_f = CompositeGradientFig4(idxFemale);

X = [X_m; X_f];
Y = [Y_m; Y_f];
group = [ones(length(X_m),1); 2*ones(length(X_f),1)];

[stats, levels] = aoctool(X, Y, group);

%% ========================================================================
%% PANEL E — Prediction Plot with 95% Prediction Interval
% ========================================================================

interaction = CompositeGradientFig4 .* TraitAnxietyFig4;

x = zscore(thickness);
y = zscore(interaction);

size_data  = 1 ./ (interaction .* 50);
color_data = zscore(TraitAnxietyFig4);

mdl = fitlm(x, y);

x_pred = linspace(min(x), max(x), 100)';
[y_pred, y_ci] = predict(mdl, x_pred, 'Prediction', 'observation');

fg3=figure('Position', [0 0 510 510]);
hold on;

scatter(x, y, size_data .* (1000000/2), color_data, 'filled');
plot(x_pred, y_pred, 'k:', 'LineWidth', 4.75);
plot(x_pred, y_ci(:,1), 'k--', 'LineWidth', 2.75);
plot(x_pred, y_ci(:,2), 'k--', 'LineWidth', 2.75);

colormap("turbo");
cb = colorbar;
ylabel(cb, 'Trait anxiety (zscore)');

xlabel({'Cortical thickness:'; ...
         'Left (PreCG + Paracentral Lobule) zscore'}, 'FontSize', 14);
ylabel({'Interaction: Att Gradient × Trait Anxiety'; ...
         '(zscore)'}, 'FontSize', 14);

set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'FontSize', 14, ...
    'XTick', -3:1:3);

xlim([-3 3]);
ylim([-3 3]);
axis square;
hold off;
exportgraphics(fg3, 'Figure4E.tif', 'Resolution', 300);
%% ========================================================================
%% PANEL F — Permutation Testing Histogram
% ========================================================================

%load('corticalthicknessInteractionWorks_31July2025.mat');

fg4=figure('Position', [0 0 350 400]);
hold on;

histogram(r_null, 'Normalization', 'pdf');
xline(empiricalCV_r, 'r:', 'LineWidth', 3);

text(empiricalCV_r - 0.1, 1, ...
    ['Empirical CV r = ' num2str(empiricalCV_r)], ...
    'FontSize', 12, ...
    'Rotation', 90);

xlabel('Null distribution of Pearson''s r');
ylabel('Density');

set(gca, 'LineWidth', 1.5, 'TickDir', 'out', ...
    'FontSize', 14, ...
    'XTick', -0.6:0.2:0.4, ...
    'YTick', 0:0.5:2.5);

hold off;
exportgraphics(fg4, 'Figure4F.tif', 'Resolution', 300);
