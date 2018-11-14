% This function calculates Binomial coefficients for choosing 0:n out of n,
% in log scale. It is useful for precalculating the coefficients and
% providing them to "Pej_Binom_test_fast.m" or "Pej_pdf_binom_fast.m" when they are redone on the same N,
% over and over again.

% Pejman, 2018, Scripps
function log_z = Pej_log_BinCoeffs(n)
% get log of all binomial coefficients for choosing 0:n out of n
log_z = nan(1,n+1);
log_z([1 n+1])=0; % the first and the last element are 0

Gammln = gammaln(1:n+1); % precalculate gamma function

m1 = floor(n/2); % this is the middle.  binomial coeffs are symmetric so we calculate only half
m2 =  ceil(n/2); % this is the "other" middle.
x  = 1:m1;
xc = n-x;
log_z(2:m1+1) = (Gammln(n+1)-Gammln(x+1)-Gammln(xc+1));

log_z(m2+1:end-1)=log_z(m1+1:-1:2); % fill in the other side
end