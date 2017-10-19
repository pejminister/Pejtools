function M = Pej_Median_withNaNs(X,d)
if length(size(X))>2
    error('I was too lazy to write this for more than 2-D')
end
if isempty(X)
    M=[];
end

if nargin>1
    if d==1
        for i = size(X,2):-1:1
            M(1,i) = Median_withNaNs_vector(X(:,i));
        end
    elseif d==2
        for i = size(X,1):-1:1
            M(i,1) = Median_withNaNs_vector(X(i,:));
        end
    else
        error('d is too big!')
    end
    
    
elseif isvector(X)
    M = Median_withNaNs_vector(X);
else
    % Get the median over the first dimension
    for i = size(X,2):-1:1
        M(1,i) = Median_withNaNs_vector(X(:,i));
    end
end
end


function M = Median_withNaNs_vector(X)
M = median(X(~isnan(X)));
end