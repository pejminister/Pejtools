% this script gives the inverse of a continous function using numerical
% optimization for a set of given points:
% y0 = fx(x0)
% => x0 = Pej_Finverse(y0, fx)

% y0 is a bunch of numbers, scalar, matrix vector whatever.
% fx is ONE function of x  like : @(x)x2+1
% x_range is [xmin, xmax] is the domain

% Note: I don't make sure that the fx is monotonic, so know that! So don't
% use this on weird functions, only smooths and easy stuff

% Pejman, New Yok, Marck 2017
% pejman.m@gmail.com

function x0 = Pej_Finverse(y0, fx, x_range)
x0 = nan(size(y0));
if nargin<3
    x_range = [-inf, +inf];
end

N = 1000;
%% quick interpolation answer
x_range(1) = max(-realmax('double'), x_range(1));
x_range(2) = min(+realmax('double'), x_range(2));

if x_range(1)>x_range(2)
    error('range(1) should be smaller than x_range(2)!')
end

if x_range(1)>0
    logdom = logspace(log10(x_range(1)), log10(x_range(2)), N);
elseif x_range(1)==0
    logdom = logspace(log10(realmin('double')), log10(x_range(2)), N);
elseif x_range(2)>0
    logdom = [-logspace(log10(-x_range(1)), log10(realmin('double')), N/2) 0 logspace(log10(realmin('double')), log10(x_range(2)), N/2)];
elseif x_range(2)==0
    logdom = -logspace(log10(-x_range(1)), log10(realmin('double')), N);
else % x_range(2)<0
    logdom = -logspace(log10(-x_range(1)), log10(-x_range(2)), N);
end

lindom = linspace(x_range(1), x_range(2), N );
xx = unique([logdom lindom]);
yy = fx(xx);
F = isfinite(yy) & imag(yy)==0;
xx = xx(F);
yy = yy(F);

for i = 1:length(y0(:))
    % do the fine tuning
    x0(i) = Pej_finverse_numerical_4one(y0(i), fx, xx,yy);
end
end


function x0 = Pej_finverse_numerical_4one(y0, fx, xx,yy)

x0_init = interp1(yy, xx, y0);
x0 = fzero( @(x)(fx(x)-y0), x0_init);
end