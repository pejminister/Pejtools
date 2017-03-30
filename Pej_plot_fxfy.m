% This function wraps matlabs plot function except that it trnasforms data
% points and axis labels using th two provided functions. This is similar
% to using gca.XScale but more general.

% For example:
% Pej_plot_fxfy(x,y,@(x)log(x), @(y)y, '.', linewidth', 2)
% would be equalt to:
% semilogx(x,y, '.', linewidth', 2);

% Example 2: plot is square root space
% Pej_plot_fxfy(x,y,@(x)sqrt(x), @(y)sqrt(y), '.')

% ----------------------------------
% Pejman, New York. March 30, 2017
% pejman.m@gmail.com
% ----------------------------------

function Pej_plot_fxfy(x,y,fx, fy, varargin)

plot(fx(x), fy(y), varargin{:});

XtickLabels = Pej_Finverse(get(gca, 'XTick'), fx);
set(gca, 'XTickLabel', XtickLabels);

YtickLabels = Pej_Finverse(get(gca, 'YTick'), fy);
set(gca, 'YTickLabel', YtickLabels);

end