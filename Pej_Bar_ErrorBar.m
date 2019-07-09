% this function plots a bar plot with error bars on it.
function h = Pej_Bar_ErrorBar(x, xl, xu, Labels, varargin)
L = length(x);
h = bar(1:L, x, 'k', varargin{:});
hold on

BV = get(h,'BaseValue');
for i =L:-1:1
if  x(i)>BV
    errorbar(i, x(i), x(i)-xl(i)  , xu(i)+nan, '.', 'linewidth', 1.5, 'color', 'w', 'Marker', 'none')
    errorbar(i, x(i), xl(i)+nan, xu(i)-x(i)  , '.', 'linewidth', 1.5, 'color', 'k', 'Marker', 'none')
elseif x(i)<BV
   errorbar(i, x(i), x(i)-xl(i)  , xu(i)+nan, '.', 'linewidth', 1.5, 'color', 'k', 'Marker', 'none')
    errorbar(i, x(i), xl(i)+nan, xu(i)-x(i)  , '.', 'linewidth', 1.5, 'color', 'w', 'Marker', 'none')
else
    errorbar(i, x(i), x(i)-xl(i)  , xu(i)+nan, '.', 'linewidth', 1.5, 'color', 'k', 'Marker', 'none')
    errorbar(i, x(i), xl(i)+nan, xu(i)-x(i)  , '.', 'linewidth', 1.5, 'color', 'k', 'Marker', 'none')
end
end
if nargin >3
    set(gca, 'XTick', 1:L, 'XTicklabel', Labels, 'XTickLabelRotation',90)
end
xlim([.4 L+.6])
set(gcf, 'position', [1 1 205  50*L]);
end