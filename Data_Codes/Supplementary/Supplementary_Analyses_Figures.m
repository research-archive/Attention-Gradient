%% =========================================================================
%% Analysis and Figures - Supplementary Information – EJN (Revision 2)
% -------------------------------------------------------------------------
%
% This script reproduces all supplementary analyses and figures.
% Required file in working directory: pearson_ci
% Toolbox: Statistics and Machine Learning Toolbox ver 12.5 and above. 
%
% Required variables inside .mat files: N = 60 
%       Neutral            [N x 3]
%       Fear               [N x 3]
%       Scrambled          [N x 3]
%       CompositeGradient  [N x 1]
%       TraitAnxiety       [N x 1]
%       StateAnxiety       [N x 1]
%       GradientNorm_First3 [N x 1]
%       GradientNorm_Last3 [N x 1]
%  % Date: 22 June 2026 
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
% Correlation between State and Trait Anxiety 
% Correlation between State Anxiety and Composite Gradient 
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


%=========================================================================
% Correlation between State and Trait Anxiety 
[r_TS,p_TS] = corr(TraitAnxiety,StateAnxiety,'Type','Pearson', 'Rows','complete');

%=========================================================================
% Correlation between State Anxiety and Composite Gradient 
[r_StateGrad,p_StateGrad] = corr(StateAnxiety,CompositeGradient,'Type','Pearson','Rows','complete');

%=========================================================================
% 95% CIs
N = sum(~isnan(TraitAnxiety) & ~isnan(StateAnxiety));
SEz = 1/sqrt(N-3);

z = atanh(r_TS);
CI_TS = tanh([z-1.96*SEz z+1.96*SEz]);

N2 = sum(~isnan(StateAnxiety) & ~isnan(CompositeGradient));
SEz2 = 1/sqrt(N2-3);

z = atanh(r_StateGrad);
CI_StateGrad = tanh([z-1.96*SEz2 z+1.96*SEz2]);

% Display

fprintf('\nTrait vs State Anxiety\n');
fprintf('r = %.3f, p = %.4f\n',r_TS,p_TS);
fprintf('95%% CI = [%.3f %.3f]\n',CI_TS(1),CI_TS(2));

fprintf('\nState Anxiety vs Composite Gradient\n');
fprintf('r = %.3f, p = %.4f\n',r_StateGrad,p_StateGrad);
fprintf('95%% CI = [%.3f %.3f]\n',CI_StateGrad(1),CI_StateGrad(2));


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

%% ========================================================================
%% Supplementary Figure 2:
% Distribution of signed Attention Gradients (Near – Far)
% =========================================================================
GradFear = Fear(:,1) - Fear(:,3);
GradNeutral = Neutral(:,1) - Neutral(:,3);
GradScrambled = Scrambled(:,1) - Scrambled(:,3);


for nParticipants = 1:60
signedGrad(nParticipants,1) = mean([GradFear(nParticipants,1), GradNeutral(nParticipants,1), GradScrambled(nParticipants,1)],2);
end 

% Fraction of positive and negative values 
prc_positive_values = ((length(signedGrad(signedGrad > 0)))/nParticipants)*100
prc_negative_values = (length(signedGrad(signedGrad < 0))/nParticipants)*100

fig = figure
histogram(signedGrad,10)
xlabel('Signed Gradient')
ylabel('Frequency')
title('Distribution of Signed Attention Gradients')
exportgraphics(fig,'Figure2Supplementary.tif','Resolution',300);


%% =========================================================================
%% Supplementary Figure 3:
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
exportgraphics(fig, 'Figure3Supplementary.tif', 'Resolution',300);

%% =========================================================================
%% Supplementary Figure 4:
%  Relationship of Cortical thickness with Attention Gradient 
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
exportgraphics(fig, 'Figure4Supplementary.tif', 'Resolution',300);


%% =========================================================================
%% Supplementary Analysis 3:
%  3A) Split-half reliability of the Euclidean CompositeGradient metric followed by Spearman-Brown correction. 
%  3B) Correlation between the Euclidean  norms of the CompositeGradient and 3-point
%  linear SlopeGradient + Scatter plot for convergent validity 
%   
% =========================================================================


%% Clear workspace
clear; close all; clc;

%% Load Data
dataFile = 'FigureSupplementary_D.mat';
load(dataFile);

%% ============================================================
% SPLIT-HALF RELIABILITY 
% ============================================================

% z-scored values of GradientNorm_First3 and GradientNorm_Last3
validIdx = ...
    ~isnan(GradientNorm_First3) & ...
    ~isnan(GradientNorm_Last3);

FirstHalf = GradientNorm_First3(validIdx);

SecondHalf = GradientNorm_Last3(validIdx);

N = length(FirstHalf);

[r_split,p_split] = corr( ...
    FirstHalf,...
    SecondHalf,...
    'Type','Pearson');

df = N - 2;

%% ============================================================
% 95% Confidence Interval
%
% Fisher z transform
% ============================================================

SEz = 1/sqrt(N-3);

z = atanh(r_split);

CI_L = tanh(z - 1.96*SEz);

CI_U = tanh(z + 1.96*SEz);

%% ============================================================
% Spearman-Brown Reliability
% ============================================================

SB = (2*r_split)/(1+r_split);

%% ============================================================
% Display Results
% ============================================================

fprintf('\n');
fprintf('=============================================\n');
fprintf('SPLIT-HALF RELIABILITY OF GRADIENT METRIC\n');
fprintf('=============================================\n');

fprintf('Pearson r(%d) = %.4f\n', ...
    df, r_split);

fprintf('95%% CI = [%.4f, %.4f]\n', ...
    CI_L, CI_U);

fprintf('p = %.6f\n', ...
    p_split);

fprintf('\nSpearman-Brown Reliability = %.4f\n', ...
    SB);

%% =======================================================================
% ============Analysis 3B + Suppementary Figure 5=========================

%% ==============================================================
% Define eccentricities
% ===============================================================

clc;
ecc = [1.5 3 6];

nParticipants = size(Neutral,1);

%% ==============================================================
% Preallocate slope variables
% ==============================================================

BetaNeutral   = zeros(nParticipants,1);
BetaFear      = zeros(nParticipants,1);
BetaScrambled = zeros(nParticipants,1);


%% ==============================================================
% R-squared values from linear fits
% ==============================================================

R2_Neutral   = zeros(nParticipants,1);
R2_Fear      = zeros(nParticipants,1);
R2_Scrambled = zeros(nParticipants,1);



%% ==============================================================
% Compute participant-specific slopes
%
% Model:
% nIES = b0 + b1(Eccentricity)
%
% Store slope coefficient (b1)
% ==============================================================

for p = 1:nParticipants

    % ---------------- Neutral ----------------
    coeffs = polyfit(ecc, Neutral(p,:), 1);

    BetaNeutral(p) = coeffs(1);

    yhat = polyval(coeffs,ecc);

    SSres = sum((Neutral(p,:) - yhat).^2);

    SStot = sum((Neutral(p,:) - mean(Neutral(p,:))).^2);

    R2_Neutral(p) = 1 - SSres/SStot;

    % ---------------- Fear -------------------
    coeffs = polyfit(ecc, Fear(p,:), 1);

    BetaFear(p) = coeffs(1);

    yhat = polyval(coeffs,ecc);

    SSres = sum((Fear(p,:) - yhat).^2);

    SStot = sum((Fear(p,:) - mean(Fear(p,:))).^2);

    R2_Fear(p) = 1 - SSres/SStot;

    % ---------------- Scrambled --------------
    coeffs = polyfit(ecc, Scrambled(p,:), 1);

    BetaScrambled(p) = coeffs(1);

    yhat = polyval(coeffs,ecc);

    SSres = sum((Scrambled(p,:) - yhat).^2);

    SStot = sum((Scrambled(p,:) - mean(Scrambled(p,:))).^2);

    R2_Scrambled(p) = 1 - SSres/SStot;

end

%% ==============================================================
% Compute linear 3 point slope-based Euclidean norm
%
% Analogous to original Near minus Far slope-based CompositeGradient
% ==============================================================

SlopeNorm = sqrt( ...
      BetaNeutral.^2 ...
    + BetaFear.^2 ...
    + BetaScrambled.^2 );


%% ==============================================================
% Participant-level mean R² across emotions
% ==============================================================

MeanR2 = mean( ...
    [R2_Neutral R2_Fear R2_Scrambled], ...
    2);

%% ==============================================================
% Recompute emotion-specific Near-Far gradients
%
% Gradient = Near - Far
%
% Near = column 1 (1.5 degrees)
% Far  = column 3 (6 degrees)
% ==============================================================

GradientNeutral = Neutral(:,1) - Neutral(:,3);

GradientFear = Fear(:,1) - Fear(:,3);

GradientScrambled = Scrambled(:,1) - Scrambled(:,3);

%% ==============================================================
% PRIMARY ANALYSIS
%
% CompositeGradient vs SlopeNorm
% ==============================================================

[r_Primary,p_Primary] = corr( ...
    CompositeGradient, ...
    SlopeNorm, ...
    'Type','Pearson');

%% ==============================================================
% SECONDARY ANALYSES
%
% Emotion-specific gradient vs slope
% ==============================================================

[r_Neutral,p_Neutral] = corr( ...
    GradientNeutral, ...
    BetaNeutral, ...
    'Type','Pearson');

[r_Fear,p_Fear] = corr( ...
    GradientFear, ...
    BetaFear, ...
    'Type','Pearson');

[r_Scrambled,p_Scrambled] = corr( ...
    GradientScrambled, ...
    BetaScrambled, ...
    'Type','Pearson');

%% ==============================================================
% Degrees of freedom
% ==============================================================

N  = nParticipants;
df = N - 2;

%% ==============================================================
% 95% Confidence Intervals for correlations
%
% Fisher z transformation
% ==============================================================

SEz = 1/sqrt(N-3);

% ---------------- Primary ----------------

z = atanh(r_Primary);

CI_Primary_L = tanh(z - 1.96*SEz);
CI_Primary_U = tanh(z + 1.96*SEz);

% ---------------- Neutral ----------------

z = atanh(r_Neutral);

CI_Neutral_L = tanh(z - 1.96*SEz);
CI_Neutral_U = tanh(z + 1.96*SEz);

% ---------------- Fear ----------------

z = atanh(r_Fear);

CI_Fear_L = tanh(z - 1.96*SEz);
CI_Fear_U = tanh(z + 1.96*SEz);

% ---------------- Scrambled ----------------

z = atanh(r_Scrambled);

CI_Scrambled_L = tanh(z - 1.96*SEz);
CI_Scrambled_U = tanh(z + 1.96*SEz);



%% ==============================================================
% Summary statistics for R²
% ==============================================================

Mean_R2   = mean(MeanR2);

Median_R2 = median(MeanR2);

SD_R2     = std(MeanR2);

SE_R2     = SD_R2/sqrt(nParticipants);

CI_R2_L = Mean_R2 - 1.96*SE_R2;

CI_R2_U = Mean_R2 + 1.96*SE_R2;
%% ==============================================================
% Display results
% ==============================================================

fprintf('\n');
fprintf('=====================================================\n');
fprintf('PRIMARY ANALYSIS\n');
fprintf('CompositeGradient vs SlopeNorm\n');
fprintf('-----------------------------------------------------\n');
fprintf('r(%d) = %.4f\n', df, r_Primary);
fprintf('95%% CI = [%.4f, %.4f]\n', ...
    CI_Primary_L, CI_Primary_U);
fprintf('p = %.6f\n', p_Primary);
fprintf('=====================================================\n');

fprintf('\n');
fprintf('SECONDARY ANALYSES\n');
fprintf('=====================================================\n');

fprintf('\nNeutral Gradient vs Neutral Slope\n');
fprintf('r(%d) = %.4f\n', df, r_Neutral);
fprintf('95%% CI = [%.4f, %.4f]\n', ...
    CI_Neutral_L, CI_Neutral_U);
fprintf('p = %.6f\n', p_Neutral);

fprintf('\nFear Gradient vs Fear Slope\n');
fprintf('r(%d) = %.4f\n', df, r_Fear);
fprintf('95%% CI = [%.4f, %.4f]\n', ...
    CI_Fear_L, CI_Fear_U);
fprintf('p = %.6f\n', p_Fear);

fprintf('\nScrambled Gradient vs Scrambled Slope\n');
fprintf('r(%d) = %.4f\n', df, r_Scrambled);
fprintf('95%% CI = [%.4f, %.4f]\n', ...
    CI_Scrambled_L, CI_Scrambled_U);
fprintf('p = %.6f\n', p_Scrambled);

fprintf('\n');
fprintf('=====================================================\n');
fprintf('LINEAR FIT QUALITY\n');
fprintf('=====================================================\n');

fprintf('Mean R^2   = %.4f\n',Mean_R2);

fprintf('Median R^2 = %.4f\n',Median_R2);

fprintf('95%% CI     = [%.4f, %.4f]\n', ...
    CI_R2_L, CI_R2_U);

%% Scatter plot - Supplementary Figure 5 

%% Figure
%% Z-score variables
x = zscore(CompositeGradient);
y = zscore(SlopeNorm);

%% Correlation
[r,p] = corr(x,y,'Type','Pearson');

%% Linear fit
mdl = fitlm(x,y);

xfit = linspace(min(x),max(x),100)';
[yfit,yCI] = predict(mdl,xfit);

%% Figure
fig = figure('Color','w');
hold on

scatter(x,y,60,'filled','MarkerFaceAlpha',0.8)

% 95% CI band
fill([xfit;flipud(xfit)],...
     [yCI(:,1);flipud(yCI(:,2))],...
     [0.65 0.65 0.65],...
     'EdgeColor','none',...
     'FaceAlpha',0.5);

plot(xfit,yfit,'k-','LineWidth',2);

%% Equal axis scaling
lims = [floor(min([x;y]))-0.5  ceil(max([x;y]))+0.5];
xlim(lims)
ylim(lims)

axis square
pbaspect([1 1 1])

xlabel('Composite Gradient (z-score)')
ylabel('SlopeNorm (z-score)')

box on
set(gca,'FontSize',12,'LineWidth',1.2)

%% Statistics text
txt = sprintf(['r(58) = 0.971\n' ...
               '95%% CI [0.951, 0.983]\n' ...
               '\\itp\\rm < .001']);

text(0.05,0.95,txt,...
    'Units','normalized',...
    'VerticalAlignment','top',...
    'FontSize',11,...
    'FontWeight','bold',...
    'Interpreter','tex');

hold off

% Export Figure
exportgraphics(fig, 'Figure5Supplementary.tif', 'Resolution',300);
