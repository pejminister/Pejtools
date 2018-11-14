function x = Pej_BetaBinomial_rnd(N, p, vScale)
if size(N,2)>1
    error('N is assumed to be a nx1 column')
end

vScale = exp(vScale);


if numel(p)==1
    p = ones(size(N))*p;
end
if numel(vScale)==1
    vScale = ones(size(N))*vScale;
end



a = vScale .* p;
b = a .* (1-p)./p;

q = betarnd(a,b);
x = binornd(N,q);
end
