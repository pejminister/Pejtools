% https://en.wikipedia.org/wiki/Logit-normal_distribution
function p = Pej_LogitNormal_pdf(x, m, s)
t = log(x./(1-x));
p = normpdf(t,m, s) ./ (x.*(1-x));

end