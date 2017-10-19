% This plots two fileds of a structure vs each other. If these's a Filter
% provided that will be used to exclude rows

function Pej_Struct_PlotFields(D,Filed1,Field2,Filter, Square, Title)
if nargin<4 || isempty(Filter)
    Filter = true(size(D.(Filed1),1),1);
end

if nargin<5
    Square = false;
end

if nargin<6
    Title='';
end

X=D.(Filed1);
Y=D.(Field2);

figure
plot(X(Filter), Y(Filter), 's', 'markersize', 3, 'markerfacecolor', [.7 0.2 0],'color', [.7 0.2 0]);
[R1,P1]=corr(X(Filter), Y(Filter));
[R2,P2]=corr(X(Filter), Y(Filter),'type', 'spearman');
tstring =  {
    sprintf('%s', strrep(Title, '_', ' '))
    sprintf('%d points, R^2:%.2f(10^{%.1f}), R^2(r):%.2f(10^{%.1f})', length(X(Filter)), R1,log10(P1), R2,log10(P2))
    };

title(tstring)
xlabel(strrep(Filed1, '_', ' '))
ylabel(strrep(Field2, '_', ' '))

set(gcf, 'position', [1 1 430         340]*.75);
if Square
    adddiag
    set(gcf, 'position', [1 1 205         200]);
    Tmpm = get(gca, 'position');
    set(gca, 'position', Tmpm([1 1 3 3]));
end

% plot regression line
Rc = regress(Y(Filter), [X(Filter) ones(size(X(Filter)))]);
% Rc = robustfit(X(Filter),Y(Filter));
hold on
plot(xlim, [xlim' ones(2,1)]*Rc, 'b', 'linewidth', 1.2)
set(gcf, 'position', [ 1   117   205   200]);

Pej_SavePlot(gcf, ['Figures/' Filed1 '_vs_' Field2 '_' Title])
end


function adddiag
M = max([ylim xlim]);
m = min([ylim xlim]);
xlim([m M])
ylim([m M])
hold on
plot([m M], [m M], 'k');
end

