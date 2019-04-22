function [RHO,PVAL] = Pej_Corr_withNaNs(X,Y, varargin)
if ~isvector(X) || ~isvector(Y)
    error('I was too lazy to write this for more than 1-D')
end
if isempty(X)
    RHO=[];
    PVAL=[];
end
F = ~isnan(X+Y);

if ~any(F)
    RHO = nan(2,2);
    PVAL= nan;
    return
end

if nargin > 2
    [RHO,PVAL] = corr(X(F), Y(F), varargin{:});
else
    [RHO,PVAL] = corr(X(F), Y(F));
end
end