%this makes a scatter plot and lowess on it
function Pej_Scatter_LOWESS(X, Y, RobustFit)
if nargin<3
    RobustFit=false;
end
F = ~isfinite(X+Y);
X(F)=[];
Y(F)=[];

[X, t] = sort(X);
Y = Y(t);

plot(X, Y, '.', 'Color', Pej_Color('blue'));
z = malowess(X, Y, 'Span', .2, 'Robust', RobustFit);

hold on
plot(X, z, '.-', 'Color', Pej_Color('red'));
end