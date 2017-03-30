function px = Pej_BetaBinomial_test(x, xc, p, vScale)
%vScale = exp(vScale);
%a = vScale * p;
%b = a * (1-p)/p;
if numel(p)==1
    p = ones(size(x))*p;
end
if numel(vScale)==1
    vScale = ones(size(x))*vScale;
end
px = zeros(size(x));
for k = 1:numel(x)
    px(k) = Pej_BetaBinomial_cdf_single(x(k), xc(k), p(k), vScale(k));
end
end

function px = Pej_BetaBinomial_cdf_single(x, xc, p, vScale)
% We assume X < Xc always
if x>xc
    t=xc;
    xc=x;
    x=t;
    p=1-p;
end

px = 0;
% px2 = 0;

vScale2 = exp(vScale);
a = vScale2 .* p;
b = a .* (1-p)./p;
t1 = gammaln(x+xc + 1);
t2 = betaln(a,b);
for i = 0:x
    txc = xc+x-i;
%         px = Pej_BetaBinomial(i, xc+x-i, p, vScale) +px;
    px = exp(t1-gammaln(i + 1)-gammaln(txc + 1)+...
        betaln((a + i),(b + txc))-t2) + px;
end
if px>.5
    px = 1-px;
end
px = px*2;% Two tailed


if px>1
    1
end
end