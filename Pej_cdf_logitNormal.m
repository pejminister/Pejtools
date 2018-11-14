% this is CDF for logit-normal distribution
% Pejman

function cx = Pej_cdf_logitNormal(x,mu,Std)
cx = nan(size(x));
F = x<1 & x>0;
cx(x<=0) = 0;
cx(x>=1) = 1;

cx(F) = .5*(1+erf((Pej_Transform_logit(x(F))-mu)/(Std*sqrt(2))));
end