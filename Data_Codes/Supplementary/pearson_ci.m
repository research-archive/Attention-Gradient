function [r, p, r_CI] = pearson_ci(x, y, alpha)
% PEARSON_CI  Computes Pearson's r, p-value, and confidence interval
%
%   [r, p, r_CI] = pearson_ci(x, y)
%   [r, p, r_CI] = pearson_ci(x, y, alpha)
%
% Inputs:
%   x, y   - numeric vectors of equal length
%   alpha  - significance level (default = 0.05 for 95% CI)
%
% Outputs:
%   r      - Pearson correlation coefficient
%   p      - two-tailed p-value
%   r_CI   - confidence interval for r [lower upper]
%
% Method:
%   Pearson correlation with CI computed using Fisher r-to-z transform.

    % Default alpha
    if nargin < 3
        alpha = 0.05;
    end

    % Remove NaNs pairwise
    valid = ~isnan(x) & ~isnan(y);
    x = x(valid);
    y = y(valid);

    % Ensure column vectors
    x = x(:);
    y = y(:);

    % Compute correlation
    [r, p] = corr(x, y, 'Type', 'Pearson');

    % Sample size
    n = length(x);

    % Fisher r-to-z transformation
    z = atanh(r);
    SE = 1 / sqrt(n - 3);
    z_crit = norminv(1 - alpha/2);

    z_CI = z + [-1 1] * z_crit * SE;
    r_CI = tanh(z_CI);

end