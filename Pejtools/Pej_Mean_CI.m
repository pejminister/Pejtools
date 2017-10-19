% This finction takes a vector X and calculates the 95CI for the mean estimate
function CI95 = Pej_Mean_CI(X)
if ~isvector(X)
    error('I was too lazy to write this for matrices, use only vector inputs : )')
end
 
if any(isnan(X))
    warning('Yo there were NaNs in the data, I discard them from the analysis')
end
X(isnan(X))=[]; % discard nans

SEM = std(X)/sqrt(length(X));               % Standard Error
ts = tinv([0.025  0.975],length(X)-1);      % T-Score
CI95 = mean(X) + ts*SEM; 
end