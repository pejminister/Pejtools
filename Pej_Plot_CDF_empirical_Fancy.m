function Pej_Plot_CDF_empirical_Fancy(X, N, HistBins,Color, NoScale)
if nargin<2 || isempty(N)
    N = 100;
end
if nargin<3 || isempty(HistBins)
    HistBins = 25;
end
if nargin<4 || isempty(Color)
    Color = [1 1 1] * 0;
end

if nargin<5 || isempty(NoScale)
    NoScale = false;
end

figure
hist(X,HistBins);
Pej_Polish_Bar

yyaxis right 
Pej_Plot_CDF_empirical(X, N, Color)
% ylim([0 YL./length(f)*100])
ylabel('Relative freq. (%)')

yyaxis left 
hist(X,HistBins);
Pej_Polish_Bar

% xlim([-.125 10.125])
% YL = max(ylim); 
% ylim([0, YL])
if NoScale==false
ylim([0, length(X)])
else
    % nothing 
end

ylabel('Count')
xlabel('X')
% set(gca, 'XTick', 0:2:10)
% set(gca, 'XTickLabel', 0:2:10)
set(gcf, 'Position', [200 200 230 200])
end