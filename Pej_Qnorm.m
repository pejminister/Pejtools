% this function Quantile normalizes the columns but keeps the median and
% MAD
% pej
function Y = Pej_Qnorm(X)
if isvector(X)
    if size(X,1)==1
        beep
        disp('Input is a row vector. This function works ONLY on columns')
    end
end
Y = nan(size(X));

for i=1:size(X,2)
    Med = median(X(:,i));
    Mad = mad(X(:,i), 1);
    
    t = tiedrank(X(:,i));
    t = (t+1)/(max(t)+2);
    Y(:,i) = norminv(t);
    
    %% salvage the scale
    Y(:,i) = Y(:,i) - median(Y(:,i));
    Y(:,i) = Y(:,i) / mad(Y(:,i),1);
    Y(:,i) = Y(:,i) * Mad;
    Y(:,i) = Y(:,i) + Med;
end
end