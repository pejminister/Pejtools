% this function Quantile normalizes the columns to standard normal.
 % pej
function Y = Pej_Qnorm(X)
Y = nan(size(X));

for i=1:size(X,2)
t = tiedrank(X(:,i));
t = (t+1)/(max(t)+2);
Y(:,i) = norminv(t);
end
end