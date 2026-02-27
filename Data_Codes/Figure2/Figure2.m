%% Figure 2 – EJN (Revision 1)
% -------------------------------------------------------------------------
% This script generates Figure 2 for the EJN manuscript.
% It plots participant-level nIES values across eccentricities 
% (1.5°, 3°, 6°) for Neutral, Fear, and Scrambled conditions,
% and overlays the group mean for each condition.
%
% Input:
%   Figure2data.mat (must contain variables:
%       Neutral    [60 x 3]
%       Fear       [60 x 3]
%       Scrambled  [60 x 3]
%
% Output:
%   Figure2.tif (300 dpi)
%
% Date: 24 Feb 2026
% -------------------------------------------------------------------------

%% Clear workspace 
clear; close all; clc;

%% Load data (expects file in current working directory)
dataFile = 'Figure2data.mat';
load(dataFile);   % Loads: Neutral, Fear, Scrambled

%% Figure parameters
eccentricityX = [2 4 6];                    % X-axis positions
eccentricityLabels = {'1.5','3','6'};       % Degree labels
xLimits = [1 7];
yLimits = [-1.5 1];
lineAlpha = 0.10;                          
nParticipants = size(Neutral,1);

%% Create figure
fg = figure;
set(fg, 'Position', [0 0 800 400]);

%% =========================
%  Plot participant data
%  =========================
for participant = 1:nParticipants

    % -------------------------
    % Neutral condition
    % -------------------------
    subplot(1,3,1); hold on;
    hNE = plot(eccentricityX, Neutral(participant,:), ...
        'k-o', ...
        'MarkerSize', 12, ...
        'MarkerFaceColor', 'k', ...
        'MarkerEdgeColor', [0.5 0.5 0.5], ...
        'LineWidth', 1.5);
    hNE.Color(4) = lineAlpha;

    % -------------------------
    % Fear condition
    % -------------------------
    subplot(1,3,2); hold on;
    hFE = plot(eccentricityX, Fear(participant,:), ...
        'b-^', ...
        'MarkerSize', 12, ...
        'MarkerFaceColor', 'b', ...
        'MarkerEdgeColor', [0.5 0.5 0.5], ...
        'LineWidth', 1.5);
    hFE.Color(4) = lineAlpha;

    % -------------------------
    % Scrambled condition
    % -------------------------
    subplot(1,3,3); hold on;
    hScr = plot(eccentricityX, Scrambled(participant,:), ...
        'm-s', ...
        'MarkerSize', 12, ...
        'MarkerFaceColor', 'm', ...
        'MarkerEdgeColor', [0.5 0.5 0.5], ...
        'LineWidth', 1.5);
    hScr.Color(4) = lineAlpha;

end

%% =========================
%  Overlay group means
%  =========================
NE_avg  = nanmean(Neutral, 1);
FE_avg  = nanmean(Fear, 1);
Scr_avg = nanmean(Scrambled, 1);

subplot(1,3,1);
plot(eccentricityX, NE_avg, ...
    'k-o', ...
    'MarkerSize', 16, ...
    'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', [0.5 0.5 0.5], ...
    'LineWidth', 3.75);

subplot(1,3,2);
plot(eccentricityX, FE_avg, ...
    'b-^', ...
    'MarkerSize', 16, ...
    'MarkerFaceColor', 'b', ...
    'MarkerEdgeColor', [0.5 0.5 0.5], ...
    'LineWidth', 3.75);

subplot(1,3,3);
plot(eccentricityX, Scr_avg, ...
    'm-s', ...
    'MarkerSize', 16, ...
    'MarkerFaceColor', 'm', ...
    'MarkerEdgeColor', [0.5 0.5 0.5], ...
    'LineWidth', 3.75);

%% =========================
%  Axis formatting
%  =========================
for sp = 1:3
    subplot(1,3,sp);
    ylim(yLimits);
    xlim(xLimits);
    set(gca, ...
        'Color', 'none', ...
        'XTick', eccentricityX, ...
        'XTickLabel', eccentricityLabels, ...
        'LineWidth', 1.5, ...
        'TickDir', 'out', ...
        'FontSize', 14);
    xlabel('Eccentricities (degrees)');
end

subplot(1,3,1);
ylabel('nIES (ms)');
title('NE');

subplot(1,3,2);
title('FE');

subplot(1,3,3);
title('Scr');

%% Export figure
exportgraphics(fg, 'Figure2.tif', 'Resolution', 300);

