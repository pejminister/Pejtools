% plots histogram with counts normalized to be percentage
function Fig = Pej_Hist_Percent(X, Bins, Plot_Median, Plot_95CI)
if nargin == 1
    Bins = 25;
end

if nargin < 3
    Plot_Median = true;
end

if nargin < 4
    Plot_95CI = false;
end

if isscalar(Bins)
    nBins = Bins;
    BinW  = range(X)/nBins;
    Bins = linspace(min(X),max(X), nBins+1) + BinW/2;
    Bins(end)=[];
end

Cbox = lines(7);
Cbox(2,:) = [0.8500 0.3250 0.0980];

% Fig = gcf; %figure;
hist(X, Bins);
h = findobj(gca,'Type','patch');
set(h, 'FaceColor',[0 0.4470 0.7410]*.5);
set(h, 'EdgeColor',[1 1 1]*.6);
set(h, 'linewidth',.5);


hold on
title(['Median:' num2str(Pej_Median_withNaNs(X),3)])
ylabel('Frequency (%)')
% set(Fig, 'position', [ 1  1 220 200])

tN = sum(isfinite(X));
Myperc = floor(max(ylim)/tN*20)*5;
set(gca, 'YTick', (0:5:Myperc)/100*tN);
if Myperc==0
    Myperc= floor(max(ylim)/tN*100);
    set(gca, 'YTick', (0:Myperc)/100*tN);
end
set(gca, 'YTickLabel', get(gca, 'YTick')/tN*100)
YL = ylim;
dy = range(YL)*.01;

if Plot_Median
plot([1 1]* Pej_Median_withNaNs(X), YL+dy*[2 -2], 'r', 'linewidth', 1, 'color', Cbox(2,:));
plot([  1]* Pej_Median_withNaNs(X), YL(1)+dy, 'r^', 'linewidth', 1, 'markerfacecolor', Cbox(2,:), 'markersize', 3, 'color', Cbox(2,:));
plot([  1]* Pej_Median_withNaNs(X), YL(2)-dy, 'rv', 'linewidth', 1, 'markerfacecolor', Cbox(2,:), 'markersize', 3, 'color', Cbox(2,:));
end
if Plot_95CI
    %        Q2_3 = 10.^ norminv([.025 .975],0, Noise.(NF{i}));
    Q2_3 = quantile(X, [.025 .975]);
    plot([1 1] * Q2_3(1), [0, max(ylim)*.75], 'r', 'linewidth', 1, 'color', Cbox(2,:));
    plot([1 1] * Q2_3(2), [0, max(ylim)*.75], 'r', 'linewidth', 1, 'color', Cbox(2,:));
    
    plot([1 ] * Q2_3(1), [max(ylim)*.75 ], 'rv', 'linewidth', 1, 'markerfacecolor', Cbox(2,:), 'markersize', 3, 'color', Cbox(2,:));
    plot([1 ] * Q2_3(2), [max(ylim)*.75 ], 'rv', 'linewidth', 1, 'markerfacecolor', Cbox(2,:), 'markersize', 3, 'color', Cbox(2,:));

    text(Q2_3(1), max(ylim)*.75, sprintf(' %.1f', Q2_3(1)), 'rotation', 90, 'color', Cbox(2,:));
    text(Q2_3(2), max(ylim)*.75, sprintf(' %.1f', Q2_3(2)), 'rotation', 90, 'color', Cbox(2,:)); 
end

plot(xlim, [0 0], 'k')
end