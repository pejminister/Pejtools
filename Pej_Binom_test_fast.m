% This function does Binomial test, with the option to provide Binomial
% coefficients to it. This is useful when the test is redone on the same N,
% over and over again.
% log_BinCoeffs can be precalculated using Pej_log_BinCoeffs(n) 

% Pejman, 2018, Scripps
function p_val = Pej_Binom_test_fast(X,N,p, log_BinCoeffs)
% here I assume x and N are scalars while p is a vector
p_val =  nan(size(p));
for i = 1:length(p)
        Bnp = Pej_pdf_binom_fast(N, p(i), log_BinCoeffs);
        p_F = Bnp <= (Bnp(X+1)*1.0000001); % I multiply the number by a little but to avoid numerical issues
        p_val(i) = sum(Bnp(p_F));
end
end


