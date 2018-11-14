% this is the most primitive numerical integration, and it's formatted like
% the matlb intergral function
% Fn is is a function handle
% Lb and Ub are the boundaries
% N is the number of Bins
% Pejman, 2018, Scripps
function P = Pej_integral(Fn,Lb,  Ub, N)
if nargin < 4
    N = 100;
end
rd = linspace(Lb,Ub,N+1);
rd2 = rd+(rd(2)-rd(1))/2;
rd2(end)=[];

f = feval(Fn, rd2);
P = sum(f)*(Ub-Lb)/N;
end
