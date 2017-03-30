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
P1(G) = binocdf(X(G)-1,N(G),p(G), 'upper');
P1(F) = binocdf(X(F)  ,N(F),p(F));
P1=P1*2;

P1(X==m)=1;
end


%% R does the test in a different way, here's the implementation:
% t = binopdf(X,N,M) * (1+1E-7);
% 
% Xt = X;
% P3 = nan(size(P1));
% for i = 1:length(X)
%     if X(i)==m(i)
%         P3(i)=1;
%     elseif X(i)<m(i)
%         y= sum(binopdf(ceil(m(i)):N(i),N(i), M(i))<=t(i));
%         P3(i) = binocdf(X(i), N(i), M(i)) + binocdf(N(i)-y, N(i), M(i), 'upper');
%     else
%         y= sum(binopdf(0:floor(m(i)),N(i), M(i))<=t(i));
%         P3(i) = binocdf(y-1, N(i), M(i)) + binocdf(X(i)-1, N(i), M(i), 'upper');
%     end
%     
% end