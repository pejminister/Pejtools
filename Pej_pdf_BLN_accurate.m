% This function calculate probability density function for Binomial
% Logit Normal distribution, that is a binomially distributed variable
% "x" with ratio "r" where  logit-transformed r is itself normally
% distributed with mean mu and variance v
% x + xc are the total number of trials in the binomial.

% x and xc should have the same dimentionality.
% mu and sigma are assumed to be scalars

% Pejman, Oct 2017
% pejman.m@gmail.com

function px = Pej_pdf_BLN_accurate(x, xc, mu, v)
minV = eps;%1E-3; % for variances less than this, binomial probability is reported. This is because for very small "v"s the calculations get numerically unstable.
px = nan(size(x));
if isrow(x)
    % make sure it's a column vector
    x = x';
    xc=xc';
end
if mu==inf
    px = double(xc==0);
    return
end
if mu==-inf
    px = double(x==0);
    return
end
if v <=minV
    % binomial
    px = binopdf(x, x+xc, lgist(mu));
    if ~isfinite(px)
        warning('Numerical failure!')
    end
    return
end
if length(mu)>1
    warning('this function is not written for vectors! make sure it works correctly!')
end

for i = 1:numel(x)
    px(i) = BLN_sub(x(i),xc(i), mu, v);
end
NaNflt = isnan(px);
if any(NaNflt)
    warning('Numerical failure, substituted with binomial pdf!')
    px(NaNflt) = binopdf(x(NaNflt), x(NaNflt)+xc(NaNflt), lgist(mu));
end
end

function mu = lgit(r)
% logit transformation
if any(r>1 | r<0)
    error('Input out of logit function domain of [0, 1]')
end
mu = log(r./(1-r));
end

function r = lgist(mu)
% logistic transformation
r = 1./(1+exp(-mu));
end

function px = BLN_sub(x,xc, mu, v)
lZ = (gammaln(x + xc + 1) - gammaln(x + 1) - gammaln(xc + 1));

% calculate the probability for one single point
dx = linspace(0,1,1000);
tpx= Pej_cdf_logitNormal(dx,mu,sqrt(v));

LB = max(dx(tpx<eps));
UB = min(dx(tpx>(1-eps)));
px = 1 ./ sqrt(2 * pi * v) * integral(@(r)fx(r, x, xc, mu, v, lZ), LB, UB, 'RelTol', 1E-4, 'AbsTol', 0);
end

function dp = fx(r, x, xc, mu, v,lZ)
dp = exp(log(r) .* (x - 1) + log(1 - r) .* (xc - 1) - (((lgit(r) - mu) .^ 2) ./ (2 * v)) + lZ);
end
