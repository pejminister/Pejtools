% This function does Binomial pdf, with the option to provide Binomial
% coefficients to it. This is useful when it is recalculated on the same N,
% over and over again. 
% log_BinCoeffs can be precalculated using Pej_log_BinCoeffs(n) 

% Pejman, 2018, Scripps
function Bin_p = Pej_pdf_binom_fast(n, p, log_BinCoeffs)
% this calculates the binomial pdf for 0:n, out of n
if p==0
    Bin_p = zeros(n+1,1);
    Bin_p(1) = 1;
    return
end

if p==1
    Bin_p = zeros(n+1,1);
    Bin_p(n+1) = 1;
    return
end

log_p = log(p);
log_pc= log(1-p);
x  = 0:n;
xc = n-x;
Bin_p = exp(log_p*x+log_pc*xc+log_BinCoeffs);
Bin_p = Bin_p / sum(Bin_p);
end


