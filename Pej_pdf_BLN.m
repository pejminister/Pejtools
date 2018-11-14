% this function is written to be fast. It might be inaccurate in cases
% where |mu| is large like >5 and x and xc are also on the extremes.
% Use Pej_pdf_BLN_accurate.m for more accuracy (it's sloW!!!!)
%
% This function calculate probability density function for Binomial
% Logit Normal distribution, that is a binomially distributed variable
% "x" with ratio "r" where  logit-transformed r is itself normally
% distributed with mean mu and variance v
% x + xc are the total number of trials in the binomial.

% x and xc should have the same dimentionality.
% mu and sigma are assumed to be scalars

% Pejman, Oct 2017
% pejman.m@gmail.com

function px = Pej_pdf_BLN(x, xc, mu, v)
    NN = 250; % number of points to integrate over
    minV = eps;% for variances less than this, binomial probability is reported. This is because for very small "v"s the calculations get numerically unstable.
    
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

    lZ = (gammaln(x+xc + 1)-gammaln(x + 1)-gammaln(xc + 1))-log(2*pi*v)*.5;
    
    % calculate the probability for one single point
    dx = linspace(0,1,NN);
    tpx= Pej_cdf_logitNormal(dx,mu,sqrt(v));
    
    LB = max(dx(tpx<eps));
    UB = min(dx(tpx>(1-eps)));
    
    rd = linspace(LB,UB,NN+1);
    
    rd2 = rd + (rd(2) - rd(1) ) / 2;
    rd2(end) = [];
    f = fx(rd2, x, xc, mu, v, lZ);
    px = exp(-log(NN/(UB-LB)) + log(sum(f, 2)));

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

function dp = fx(r, x, xc, mu, v, z)
    dp = exp(log(r) .* (x - 1) + log(1 - r) .* (xc - 1) - (((lgit(r) - mu) .^ 2) ./ (2 * v)) + z);
end
