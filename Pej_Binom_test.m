% this does a binomial twosided test
function P1 =Pej_Binom_test(X,N,p)
m = p.*N;
if numel(p)==1
    p = ones(size(X))*p;
end
if numel(N)==1
    N = ones(size(X))*N;
end

F = X<=m;G=~F;

P1 = nan(size(X));
if ~isempty(G)
    P1(G) = binocdf(X(G)-1,N(G),p(G), 'upper');
end
if ~isempty(F)
    P1(F) = binocdf(X(F)  ,N(F),p(F));
end
P1=P1*2;

P1(X==m)=1;
end


%% R does the test in a different way, here's the implementation:
% t = binopdf(X,N,M) * (1+1E-7);
% P1 =  nan(size(X));
% for i = 1:length(X)
%     if X(i)==m(i)
%         P1(i)=1;
%     elseif X(i)<m(i)
%         y= sum(binopdf(ceil(m(i)):N(i),N(i), M(i))<=t(i));
%         P1(i) = binocdf(X(i), N(i), M(i)) + binocdf(N(i)-y, N(i), M(i), 'upper');
%     else
%         y= sum(binopdf(0:floor(m(i)),N(i), M(i))<=t(i));
%         P1(i) = binocdf(y-1, N(i), M(i)) + binocdf(X(i)-1, N(i), M(i), 'upper');
%     end
%
% end