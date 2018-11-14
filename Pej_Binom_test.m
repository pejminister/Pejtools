% this does a binomial two-sided test
function p_val2 = Pej_Binom_test(X,N,p)
m = p.*N;
[errorcode, X, N, p] = distchck(3,X,N,p);
if errorcode > 0
    error(message('stats:binopdf:InputSizeMismatch'));
end

t = binopdf(X,N,p) * (1+1E-7);
% p_val =  nan(size(X));
p_val2 = nan(size(X));
for i = 1:length(X)
    if X(i)==m(i)
%         p_val(i)=1;
        p_val2(i)=1;
    elseif X(i)<m(i)
%         y  = sum(binopdf(ceil(m(i)):N(i),N(i), p(i))<=t(i));
%         p_val(i) = binocdf(X(i), N(i), p(i)) + binocdf(N(i)-y, N(i), p(i), 'upper');
        tp = binopdf(ceil(m(i)):N(i),N(i), p(i));
        p_othertail = sum(tp(tp<=t(i)));
        p_val2(i) =p_othertail+ binocdf(X(i), N(i), p(i));
    else
%         y= sum(binopdf(0:floor(m(i)),N(i), p(i))<=t(i));
%         p_val(i) = binocdf(y-1, N(i), p(i)) + binocdf(X(i)-1, N(i), p(i), 'upper');        
        tp = binopdf(0:floor(m(i)),N(i), p(i));
        p_othertail = sum(tp(tp<=t(i)));
        p_val2(i) =p_othertail+ binocdf(X(i)-1, N(i), p(i), 'upper');
    end

end
end