%% Figure 3 – EJN (Revision 1)
% -------------------------------------------------------------------------
% This script generates Figure 3 for the manuscript.
% It reproduces all panels associated with Figure 3, including:
%   Panel B – Correlation between right Cerebellum Lobule VI rGMV 
%               and composite attentional gradient
%   Panel D – Bilateral Cerebellum Lobule VI × Trait Anxiety interaction
%   Gender-based moderation analyses
%   Panel E – Prediction plot with 95% prediction interval
%   Panel F – Permutation testing histogram
%
% Input:
%   Fig3_1.mat  
%
%   Fig3_2.mat  
%
%   Fig3_3.mat  
%
% Output:
%   Figure 3 panels (B, D, E, F); 300 dpi 
%   Statistical values printed to MATLAB Command Window
%
% Date: 24 Feb 2026
% -------------------------------------------------------------------------

%% Clear workspace 
clearvars;
close all;
clc;

%% ========================================================================
%% PANEL B — Scatter Plot: Right Cerebellum Lobule VI
% ========================================================================

% Load data
dataFile = 'Fig3_1.mat';
load(dataFile);



% Create figure
fg1 = figure('Position', [0 0 400 400]);
hold on;

scatter( ...
    zscore(rtCb6), ...
    zscore(CompositeGradient), ...
    20*4.5, ...
    'MarkerEdgeColor', [0.5 0.5 0.5], ...
    'MarkerFaceColor', [0.2 0.2 0.2], ...
    'LineWidth', 1.25);

% Add least-squares line
ls = lsline;
ls.LineStyle = ':';
ls.LineWidth = 4.75;
ls.XData     = [-2.3 2.3];
ls.Color     = [0 0 0];

xlabel('rGMV (z score)', 'FontSize', 14);
ylabel('Att gradient (z score)', 'FontSize', 14);
title('Rt Cerebellum Lobule VI');

set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'FontSize', 14);
xlim([-2.5 2.5]);
ylim([-2.5 2.5]);
axis square;
hold off;

% Export figure
exportgraphics(fg1, 'Figure3B.tif', 'Resolution', 300);
%% ---------- Correlation statistics (reported in manuscript)

X = rtCb6;
Y = CompositeGradient;

n  = length(X);
df = n - 2;

[r, p] = corr(X, Y, 'Type', 'Pearson');

% 95% CI via Fisher transformation
z       = atanh(r);
se      = 1 / sqrt(n - 3);
z_crit  = 1.96;
z_CI    = [z - z_crit*se, z + z_crit*se];
r_CI    = tanh(z_CI);

fprintf('\n=== Right Cerebellum Lobule VI ===\n');
fprintf('r = %.4f\n', r);
fprintf('df = %d\n', df);
fprintf('95%% CI = [%.4f, %.4f]\n', r_CI(1), r_CI(2));
fprintf('p = %.4f\n\n', p);

clear X Y

%% ========================================================================
%% PANEL D — Interaction: Bilateral Cerebellum Lobule VI × Anxiety
% ========================================================================

load('Fig3_2.mat');

fg2 = figure('Position', [0 0 400 400]);
hold on;

% High anxiety
scatter( ...
    zscore(Cb6(idxHA)), ...
    zscore(CompositeGradient(idxHA)), ...
    20*4.5, ...
    'MarkerEdgeColor', [0.5 0.5 0.5], ...
    'MarkerFaceColor', [1 0 0], ...
    'LineWidth', 1.25);

% Low anxiety
scatter( ...
    zscore(Cb6(idxLA)), ...
    zscore(CompositeGradient(idxLA)), ...
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
%title('Bilat Cerebellum Lobule VI × Anxiety');

set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'FontSize', 14);
xlim([-2.5 2.5]);
ylim([-2.5 2.5]);
axis square;
hold off;

% Export figure
exportgraphics(fg2, 'Figure3D.tif', 'Resolution', 300);
%% ---------- Slopes for interaction

X1 = Cb6(idxHA); Y1 = CompositeGradient(idxHA);
X2 = Cb6(idxLA); Y2 = CompositeGradient(idxLA);

coef1 = polyfit(X1, Y1, 1);
coef2 = polyfit(X2, Y2, 1);

fprintf('Slope (High Anxiety) = %.4f\n', coef1(1));
fprintf('Intercept (High Anxiety) = %.4f\n', coef1(2));
fprintf('Slope (Low Anxiety) = %.4f\n', coef2(1));
fprintf('Intercept (Low Anxiety) = %.4f\n\n', coef2(2));

%% ========================================================================
%% Gender-Based Moderation Analysis
% ========================================================================

X_male   = Cb6(idxMale);
Y_male   = CompositeGradient(idxMale);
X_female = Cb6(idxFemale);
Y_female = CompositeGradient(idxFemale);

X = [X_male; X_female];
Y = [Y_male; Y_female];

Group = categorical([ ...
    ones(length(X_male),1); ...
    2*ones(length(X_female),1)]);

T = table(X, Y, Group);

mdl = fitlm(T, 'Y ~ X*Group');

fprintf('\n=== Gender Interaction Model ===\n');
disp(mdl);

anovaTable = anova(mdl, 'summary');
disp(anovaTable);

%% ========================================================================
%% PANEL E — Prediction Plot with 95% Prediction Interval
% ========================================================================

interaction = CompositeGradient .* TraitAnxiety;

x = zscore(Cb6);
y = zscore(interaction);

temp       = interaction .* 50;
size_data  = 1 ./ temp;
color_data = zscore(TraitAnxiety);

mdl = fitlm(x, y);

x_pred = linspace(min(x), max(x), 100)';
[y_pred, y_ci] = predict(mdl, x_pred, 'Prediction', 'observation');

fg3 = figure('Position', [0 0 510 510]);
hold on;

scatter(x, y, size_data .* (1000000/2), color_data, 'filled');
plot(x_pred, y_pred, 'k:', 'LineWidth', 4.75);
plot(x_pred, y_ci(:,1), 'k--', 'LineWidth', 2.75);
plot(x_pred, y_ci(:,2), 'k--', 'LineWidth', 2.75);

colormap("turbo");
cb = colorbar;
ylabel(cb, 'Trait anxiety (zscore)');

xlabel('rGMV: Bilateral Cerebellum Lobule VI (zscore)', 'FontSize', 14);
ylabel({'Interaction: Att Gradient × Trait Anxiety'; '(zscore)'}, 'FontSize', 14);

set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'FontSize', 14);
xlim([-3 3]);
ylim([-3 3]);
axis square;
hold off;
% Export figure
exportgraphics(fg3, 'Figure3E.tif', 'Resolution', 300);
%% ========================================================================
%% PANEL F — Permutation Testing Histogram
% ========================================================================

load('Fig3_3.mat');

fg4 = figure('Position', [0 0 350 400]);
hold on;

histogram(r_null, 'Normalization', 'pdf');
xline(empiricalCV_r, 'r:', 'LineWidth', 3);

text(empiricalCV_r - 0.1, 1, ...
    ['Empirical CV r = ' num2str(empiricalCV_r)], ...
    'FontSize', 12, ...
    'Rotation', 90);

xlabel('Null distribution of Pearson''s r');
ylabel('Density');

set(gca, ...
    'LineWidth', 1.5, ...
    'TickDir', 'out', ...
    'FontSize', 14, ...
    'XTick', -0.6:0.2:0.4, ...
    'YTick', 0:0.5:2.5);

hold off;


% Export figure
exportgraphics(fg4, 'Figure3F.tif', 'Resolution', 300);