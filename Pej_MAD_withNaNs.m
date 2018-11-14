function M = Pej_MAD_withNaNs(X,flag, d)
if length(size(X))>2
    error('I was too lazy to write this for more than 2-D')
end
if isempty(X)
    M=[];
end
if nargin<2 || isempty(flag)
    flag = 0;
end
if nargin>2
    if d==1
        for i = size(X,2):-1:1
            M(1,i) = MAD_withNaNs_vector(X(:,i),flag);
        end
    elseif d==2
        for i = size(X,1):-1:1
            M(i,1) = MAD_withNaNs_vector(X(i,:),flag);
        end
    else
        error('d is too big!')
    end
    
    
elseif isvector(X)
    M = MAD_withNaNs_vector(X,flag);
else
    % Get the mean over the first dimension
    for i = size(X,2):-1:1
        M(1,i) = MAD_withNaNs_vector(X(:,i),flag);
    end
end
end


function M = MAD_withNaNs_vector(X,flag)
M = mad(X(~isnan(X)),flag);
end