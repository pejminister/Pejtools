%this makes a series of boxplots showing the trands in X as a function of
%group Var
function Pej_Boxplot_trend(X, groupVar, nBox)
if nargin<3
    nBox=5;
end

% discurd nans
F = isnan(X+groupVar);
X(F)=[];
groupVar(F)=[];


[~, I] = sort(groupVar, 'ascend');
X = X(I);
groupVar = groupVar(I);

BinSize = (.1+length(groupVar))/nBox;
Groups = nan(size(X));
Groups(:)  = floor((1:length(groupVar))/BinSize);

Pej_BoxPlot(X, Groups)
hold on
plot([xlim], median(X)* [1 1], '--k')


