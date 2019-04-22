function [ b, sigma2_x, x_est, y_est, stats] = Pej_Deming_regression(x,y,varargin)

F = ~isnan(x+y);
x_est = nan(size(x));
y_est = nan(size(x));

[ b sigma2_x x_est(F) y_est(F) stats] = deming(x(F),y(F),varargin{:});
end