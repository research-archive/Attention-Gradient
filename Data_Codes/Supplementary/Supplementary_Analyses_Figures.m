%% =========================================================================
%% Analysis and Figures - Supplementary Information – EJN (Revision 1)
% -------------------------------------------------------------------------
%
% This script reproduces all supplementary analyses and figures.
% Required file in working directory:
% Figure2data.mat
%
% Required variables inside .mat file:
%       Neutral            [N x 3]
%       Fear               [N x 3]
%       Scrambled          [N x 3]
%       CompositeGradient  [N x 1]
%       TraitAnxiety       [N x 1]
%  % Date: 24 Feb 2026
% =========================================================================

%% Clear workspace
clear; close all; clc;

%% Load Data
dataFile = 'FigureSupplementary_A.mat';
load(dataFile);

nParticipants = size(Neutral,1);

%% =========================================================================
%% Supplementary Analysis 1:
% One-Way Repeated-Measures ANOVA Across Eccentricities
% =========================================================================

%% Compute eccentricity-wise averages (collapsed across emotion)

ECC_near = nanmean([Neutral(:,1), Fear(:,1), Scrambled(:,1)], 2);
ECC_intermediate = nanmean([Neutral(:,2), Fear(:,2), Scrambled(:,2)], 2);
ECC_far = nanmean([Neutral(:,3), Fear(:,3), Scrambled(:,3)], 2);

%% Normality (Lilliefors Test)

[H_near,P_near] = lillietest(ECC_near); % 1.5 deg eccentricity
[H_int,P_int]   = lillietest(ECC_intermediate); % 3 deg eccentricity
[H_far,P_far]   = lillietest(ECC_far); % 6 deg eccentricity

fprintf('\nNormality Tests (Lilliefors):\n');
fprintf('Near: H=%d, p=%.4f\n',H_near,P_near);
fprintf('Intermediate: H=%d, p=%.4f\n',H_int,P_int);
fprintf('Far: H=%d, p=%.4f\n',H_far,P_far);

%% Prepare Repeated-Measures Model

T = table(ECC_near, ECC_intermediate, ECC_far, ...
    'VariableNames',{'Near','Intermediate','Far'});

WithinDesign = table(categorical({'Near';'Intermediate';'Far'}), ...
    'VariableNames',{'Eccentricity'});

rm = fitrm(T,'Near-Far ~ 1','WithinDesign',WithinDesign);
ranovatbl = ranova(rm,'WithinModel','Eccentricity');

disp('Repeated-Measures ANOVA:')
disp(ranovatbl)

%% Mauchly's Test for Sphericity

mauchlyTbl = mauchly(rm);
fprintf('\nMauchly''s Test: W=%.4f, p=%.4f\n', ...
    mauchlyTbl.W, mauchlyTbl.pValue);

sphericityViolated = mauchlyTbl.pValue < 0.05;

%% Extract Statistics (Greenhouse-Geisser corrected if spericity is violated)

rowEffect = strcmp(ranovatbl.Row,'(Intercept):Eccentricity');
rowError  = strcmp(ranovatbl.Row,'Error(Eccentricity)');

Fvalue = ranovatbl.F(rowEffect);
df1    = ranovatbl.DF(rowEffect);
df2    = ranovatbl.DF(rowError);

if sphericityViolated
    pValue = ranovatbl.pValueGG(rowEffect);
else
    pValue = ranovatbl.pValue(rowEffect);
end

fprintf('\nRM-ANOVA Result:\nF(%d,%d)=%.4f, p=%.4f\n',df1,df2,Fvalue,pValue);

%% Effect size: Partial Eta Squared

SS_effect = ranovatbl.SumSq(rowEffect);
SS_error  = ranovatbl.SumSq(rowError);
eta_p2 = SS_effect/(SS_effect+SS_error);

fprintf('Partial eta squared (ηp²)=%.4f\n',eta_p2);

%% =========================================================================
%% Supplementary Analysis 2:
% Correlation Between nIES and Trait Anxiety
% =========================================================================

predictorVars = {
    Neutral(:,1); Neutral(:,2); Neutral(:,3);
    Fear(:,1); Fear(:,2); Fear(:,3);
    Scrambled(:,1); Scrambled(:,2); Scrambled(:,3);
    CompositeGradient };

predictorNames = {
    'Neutral-Near'; 'Neutral-Intermediate'; 'Neutral-Far';
    'Fear-Near'; 'Fear-Intermediate'; 'Fear-Far';
    'Scrambled-Near'; 'Scrambled-Intermediate'; 'Scrambled-Far';
    'CompositeGradient'};

numPredictors = numel(predictorVars);
storeVals = zeros(numPredictors,4);

y = TraitAnxiety;

for ii = 1:numPredictors
    [r,p,CI] = pearson_ci(predictorVars{ii}, y);
    storeVals(ii,:) = [r p CI(1) CI(2)];
end

resultsTable = table(predictorNames, ...
    storeVals(:,1), storeVals(:,2), ...
    storeVals(:,3), storeVals(:,4), ...
    'VariableNames',{'Predictor','Pearsons_r','p_value','CI_Lower','CI_Upper'});

disp('Correlation Results:')
disp(resultsTable)

writetable(resultsTable,'Correlation_Results_nIES_STAIT.csv');

%% =========================================================================
%% Supplementary Figure 1:
% Emotion-Specific Attention Gradients (Near – Far)
% =========================================================================

% Calculate emotion-specific gradients
gradient_NE  = Neutral(:,1) - Neutral(:,3);
gradient_FE  = Fear(:,1) - Fear(:,3);
gradient_Scr = Scrambled(:,1) - Scrambled(:,3);


% Create figure
fig = figure('Color','w','Position',[100 100 450 450]);
boxplot([gradient_NE,gradient_FE,gradient_Scr],'Symbol','');

set(gca,'LineWidth',1.5,...
    'FontSize',14,...
    'TickDir','out',...
    'XTickLabel',{'Neutral','Fear','Scrambled'});

ylabel('Attention gradient');
xlabel('Emotion Condition');
box off;

exportgraphics(fig,'Figure1supplementary.tif','Resolution',300);

%% =========================================================================
%% Supplementary Figure 1 Analysis:
% One-Way Repeated-Measures ANOVA Across Emotion Gradients
% =========================================================================

% Normality
[H_NE,P_NE]   = lillietest(gradient_NE);
[H_FE,P_FE]   = lillietest(gradient_FE);
[H_Scr,P_Scr] = lillietest(gradient_Scr);

fprintf('\nNormality (Gradients):\n');
fprintf('Neutral: H=%d, p=%.4f\n',H_NE,P_NE);
fprintf('Fear: H=%d, p=%.4f\n',H_FE,P_FE);
fprintf('Scrambled: H=%d, p=%.4f\n',H_Scr,P_Scr);


% Prepare table for One-Way Rep Measures ANOVA
T = table(gradient_NE,gradient_FE,gradient_Scr, ...
    'VariableNames',{'Neutral','Fear','Scrambled'});

WithinDesign = table(categorical({'Neutral';'Fear';'Scrambled'}), ...
    'VariableNames',{'Emotion'});

rm = fitrm(T,'Neutral-Scrambled ~ 1','WithinDesign',WithinDesign);
ranovatbl = ranova(rm,'WithinModel','Emotion');

disp('Gradient RM-ANOVA:')
disp(ranovatbl)

% Mauchly's test for spericity
mauchlyTbl = mauchly(rm);
fprintf('\nMauchly''s Test: W=%.4f, p=%.4f\n', ...
    mauchlyTbl.W, mauchlyTbl.pValue);

rowEffect = strcmp(ranovatbl.Row,'(Intercept):Emotion');
rowError  = strcmp(ranovatbl.Row,'Error(Emotion)');

Fvalue = ranovatbl.F(rowEffect);
df1    = ranovatbl.DF(rowEffect);
df2    = ranovatbl.DF(rowError);

if mauchlyTbl.pValue < 0.05
    pValue = ranovatbl.pValueGG(rowEffect);
else
    pValue = ranovatbl.pValue(rowEffect);
end

fprintf('\nF(%d,%d)=%.4f, p=%.4f\n',df1,df2,Fvalue,pValue);

SS_effect = ranovatbl.SumSq(rowEffect);
SS_error  = ranovatbl.SumSq(rowError);
eta_p2 = SS_effect/(SS_effect+SS_error);

fprintf('Partial eta squared (ηp²)=%.4f\n',eta_p2);

%% =========================================================================
%% Supplementary Figure 2:
%  Relationship of bilateral GMV with Attention Gradient 
% =========================================================================
%% Clear workspace
clear; close all; clc;

%% Load Data
dataFile = 'FigureSupplementary_B.mat';
load(dataFile);

% Identify gender indices
idxMale   = gender == 1;
idxFemale = gender == 2;

% Median split separately for each gender
medianMale   = median(TraitAnxiety(idxMale));
medianFemale = median(TraitAnxiety(idxFemale));

% High / Low Anxiety Groups
idxMaleHA   = idxMale   & TraitAnxiety >= medianMale;
idxMaleLA   = idxMale   & TraitAnxiety <  medianMale;

idxFemaleHA = idxFemale & TraitAnxiety >= medianFemale;
idxFemaleLA = idxFemale & TraitAnxiety <  medianFemale;

%% ---------------------------------------------------------
%% Extract Variables
%% ---------------------------------------------------------

% High Anxiety
X1 = Cb6(idxMaleHA);      Y1 = CompositeGradient(idxMaleHA);
X2 = Cb6(idxFemaleHA);    Y2 = CompositeGradient(idxFemaleHA);

% Low Anxiety
X3 = Cb6(idxMaleLA);      Y3 = CompositeGradient(idxMaleLA);
X4 = Cb6(idxFemaleLA);    Y4 = CompositeGradient(idxFemaleLA);

%% ---------------------------------------------------------
%% Create Figure
%% ---------------------------------------------------------

fig = figure('Color','w','Position',[100 100 1100 600]);

%% =========================================================
%% PANEL A – HIGH ANXIETY
%% =========================================================
subplot(1,2,1); hold on;

% Z-score for visualization only
plot(zscore(X1), zscore(Y1), 'ro', 'MarkerSize',12,'LineWidth',1.25);
plot(zscore(X2), zscore(Y2), 'r^', 'MarkerSize',12,'LineWidth',1.25);

% Regression lines
ll = lsline;
set(ll(1),'LineStyle','-','LineWidth',2.5,'XData',[-2.3 2.3]);
set(ll(2),'LineStyle',':','LineWidth',2.5,'XData',[-2.3 2.3]);

% Correlations and slopes (raw values)
[r1, ~] = corr(X1, Y1, 'Type','Pearson');
[r2, ~] = corr(X2, Y2, 'Type','Pearson');

slope1 = polyfit(X1, Y1, 1); slope1 = slope1(1);
slope2 = polyfit(X2, Y2, 1); slope2 = slope2(1);

df1 = length(X1) - 2;
df2 = length(X2) - 2;

% Statistical annotation
text(-2.4, 2.3, sprintf('Male  r(%d)=%.2f, m=%.2f', df1, r1, slope1));
text(-2.4, 2.0, sprintf('Female  r(%d)=%.2f, m=%.2f', df2, r2, slope2));

% Axis formatting
xlabel('rGMV (z score)');
ylabel('Att gradient (z score)');
xlim([-2.5 2.5]); ylim([-2.5 2.5]);
axis square;
set(gca,'LineWidth',1.5,'TickDir','out');

% Panel label
text(-3.45, 2.45, 'A', 'FontSize',24, 'FontWeight','bold');

hold off;

%% =========================================================
%% PANEL B – LOW ANXIETY
%% =========================================================
subplot(1,2,2); hold on;

plot(zscore(X3), zscore(Y3), 'bo', 'MarkerSize',12,'LineWidth',1.25);
plot(zscore(X4), zscore(Y4), 'b^', 'MarkerSize',12,'LineWidth',1.25);

ll = lsline;
set(ll(1),'LineStyle','-','LineWidth',2.5,'XData',[-2.3 2.3]);
set(ll(2),'LineStyle',':','LineWidth',2.5,'XData',[-2.3 2.3]);

[r3, ~] = corr(X3, Y3, 'Type','Pearson');
[r4, ~] = corr(X4, Y4, 'Type','Pearson');

slope3 = polyfit(X3, Y3, 1); slope3 = slope3(1);
slope4 = polyfit(X4, Y4, 1); slope4 = slope4(1);

df3 = length(X3) - 2;
df4 = length(X4) - 2;

text(-2.4, 2.3, sprintf('Male  r(%d)=%.2f, m=%.2f', df3, r3, slope3));
text(-2.4, 2.0, sprintf('Female  r(%d)=%.2f, m=%.2f', df4, r4, slope4));

xlabel('rGMV (z score)');
ylabel('Att gradient (z score)');
xlim([-2.5 2.5]); ylim([-2.5 2.5]);
axis square;
set(gca,'LineWidth',1.5,'TickDir','out');

% Panel label
text(-3.45, 2.45, 'B', 'FontSize',24, 'FontWeight','bold');

hold off;

% Export Figure
exportgraphics(fig, 'Figure2_Supplementary.tif', 'Resolution',300);

%% =========================================================================
%% Supplementary Figure 3:
%  Relationship of Cortical thickness with Attention Gradient [check not
%  working as reported !!] % note that gender and Composite gradient order is different here 
% =========================================================================
%% Clear workspace
clear; close all; clc;

% Load data 
dataFile = 'FigureSupplementary_C.mat';
load(dataFile);% Identify gender indices
idxMale   = gender == 1;
idxFemale = gender == 2;


%% Extract Variables
%% ---------------------------------------------------------

% High Anxiety
X1 = thickness(idxMale);      Y1 = thicknessCompositeGradient(idxMale);
X2 = thickness(idxFemale);    Y2 = thicknessCompositeGradient(idxFemale);

%% ---------------------------------------------------------
%% Create Figure
%% ---------------------------------------------------------

fig = figure('Color','w','Position',[100 100 400 400]);

hold on;

% Z-score for visualization only
plot(zscore(X1), zscore(Y1), 'ko', 'MarkerSize',12,'LineWidth',1.25);
plot(zscore(X2), zscore(Y2), 'ks', 'MarkerSize',12,'LineWidth',1.25);

% Regression lines
ll = lsline;
set(ll(1),'LineStyle',':','LineWidth',2.5,'XData',[-2.3 2.3]);
set(ll(2),'LineStyle','-','LineWidth',2.5,'XData',[-2.3 2.3]);

% Axis formatting
xlabel('Cortical thickness (z score)');
ylabel('Att gradient (z score)');
xlim([-2.5 2.5]); ylim([-2.5 2.5]);
axis square;
set(gca,'LineWidth',1.5,'TickDir','out');

hold off 
% Export Figure
exportgraphics(fig, 'Figure3_Supplementary.tif', 'Resolution',300);
