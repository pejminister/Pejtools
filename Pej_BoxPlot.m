function Pej_BoxPlot(X, K)
if size(X,1)==1; X = X';end

if nargin==1
L = size(X,1);
boxplot(X,'colors', [0 0.4470 0.7410]*.5 , 'symbol','')
hold on
for i = 1:size(X,2)
    tmpX= min(max(randn(L,1)/15,-.35),.35) + i;
    plot(tmpX,X(:,i), 's', 'color', [1 1 1]*.4, 'markerfacecolor', [1 1 1]*.4, 'markersize', 3)
end
else
% Use second input as a grouping variable similar to Boxplot function from
% MAtlab
[Ks,Ks_ia,Ks_ic] =  unique(K);
boxplot(X,Ks_ic, 'colors', [0 0.4470 0.7410]*.5, 'symbol','')
hold on
D = length(Ks);
for i = 1:D
%     F = K==Ks(i);
     F = Ks_ic==i;
    L = sum(F);
    if L>0
    tmpX= min(max(randn(L,1)/15,-.35),.35) + i;
    plot(tmpX,X(F), 's', 'color', [1 1 1]*.4, 'markerfacecolor', [1 1 1]*.4, 'markersize', 3)
    end
end   
boxplot(X,Ks_ic, 'colors', [ 0.8500    0.3250    0.0980], 'symbol','')    
set(gca, 'XTickLabel', Ks);
end