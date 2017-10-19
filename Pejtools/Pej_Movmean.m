% if k is even the window for nth WILL NOT include the x(n), if k is odd
% window will be symmetric around the nth value. sides are discarded.

function M = Pej_Movmean(x,k)
L = length(x);
M = nan(size(x));

if mod(k,2)==1
    % k is odd
    w = (k-1)/2;
    M(1+w)= sum(x(1:k));
   for wc =  2+w:L-w
       M(wc) = x(wc+w)-x(wc-w-1) + M(wc-1);
   end
   M = M / k;
else
    
    
end



end