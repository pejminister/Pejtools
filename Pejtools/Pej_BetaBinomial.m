function px = Pej_BetaBinomial(x, xc, p, vScale)
vScale = exp(vScale);

if numel(p)==1
    p = ones(size(x))*p;
end
if numel(vScale)==1
    vScale = ones(size(x))*vScale;
end

a = vScale .* p;
b = a .* (1-p)./p;
px = exp(gammaln(x+xc + 1)-gammaln(x + 1)-gammaln(xc + 1)+...
    betaln((a + x),(b + xc))-betaln(a,b));
end