% This is the same as corr2 MAtlab function, except that it also does rank
% corr
function c = Pej_Corr2(X, Y, Type)
if nargin<3
    Type = 'Pearson';
end
t = corr(X, Y, 'type', Type);
c = t(1);
end