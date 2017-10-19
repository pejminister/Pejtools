function Fig = Pej_Hist_Categorical(X, YScale)
if nargin < 2
    YScale = 'linear';
end
    
[T, ~, bi] = unique(X);
Fig = figure;
Mt = max(bi);

[n, xout] = hist(bi, 1:Mt);
bar(xout, n, 'barwidth', 1, 'basevalue', 1);
set(gca,'YScale',YScale)

set(gca, 'XTick', 1:Mt)
set(gca, 'XTickLabel', strrep(T, '_', ' '))
set(gca, 'XTickLabelRotation', 90)
grid on
xlim([0 Mt+1])
end