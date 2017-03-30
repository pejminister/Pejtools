function Pej_Plot_CDF_empirical(X, N)
if nargin<2
    N = 100;
end
InfSwitch = N==inf;% do not remove the end points

N = min(N, length(X));
D = linspace(0,1, N+2);

if ~InfSwitch
    D(1)=[];
    D(end)=[];
end

% if N>=100
D100 = D ;%.01:.01:.99;
% end

Qs = quantile(X,D);
Qs100 = quantile(X,D100);
plot(Qs, D*100,       'k-', 'markerfacecolor', 'k', 'markersize', 3, 'color', 'k', 'linewidth', 1); hold on
plot(Qs100, D100*100, 'sk', 'markerfacecolor', 'k', 'markersize', 3, 'color', 'k', 'linewidth', 1);
ylabel('Percentile (%)')
set(gcf, 'position', [1 1 210 200]);
end