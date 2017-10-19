function Pej_Plot_sqrtsqrt(X,Y, varargin)
plot(sqrt(X),sqrt(Y), varargin{:});

tmpTicksX = get(gca, 'XTick');
tmpTicksY = get(gca, 'YTick');
% 
% set(gca, 'XTick', tmpTicksX);
% set(gca, 'YTick', tmpTicksY);

set(gca, 'XTickLabel', tmpTicksX.^2);
set(gca, 'YTickLabel', tmpTicksY.^2);

end