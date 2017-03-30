function BIC = Pej_BIC(Fit_residuals, Nparams)
N = length(Fit_residuals);
Rstd = std(Fit_residuals);
LL = sum(log(normpdf(Fit_residuals, 0, Rstd)));

BIC = -2*LL+ Nparams*log(N);
end