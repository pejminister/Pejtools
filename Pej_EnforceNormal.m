% This function normalizes each column of X to a normal distribution using
% QQ-normalization.
% The normalized distribution will have the same sample Median and Median Abs. Dev. as the original data.

% This function will ignore NaNs in the data

% Pejman 2016
% Pejman.m@gmail.com
% -------------------------

function X = Pej_EnforceNormal(X)
if length(size(X))>2
    error('This code is written for 1D and 2D arrays only.')
end
if isvector(X)
    F = ~isnan(X);
    X(F) = Pej_EnforceNormal_v(X(F));
else
    for i = 1:size(X,2)
        F = ~isnan(X(:,i));
        X(F, i) = Pej_EnforceNormal_v(X(F, i));
    end
end

end


function Xn = Pej_EnforceNormal_v(X)
if isempty(X)
    Xn = X;
    return
end

%% break ties
[tmp_Iu,~, tmpI]= unique(X); % find ties
tmpdiff = diff(sort(X, 'ascend'));
dX = min(tmpdiff(tmpdiff>0))*1E-6;
X = X+(rand(size(X))-.5)*dX;
% Endof break ties

if sum(isfinite(X))<= (length(X)/2)
    warning('Too many infinite values, Normalization Failed! The values were not normalized.')
    Xn = X;
    return
end
M  = median(X);
Md = mad(X,1); % Median abs dev
I  = tiedrank(X);
Q  = I/(length(X)+1);
Xn = norminv(Q, 0, 1); % Make data from standard normal
Mdn =mad(Xn);
if Mdn ==0
    Mdn = 1;
end

if isnan(Md)
    Md = 0;
end
Xn = Xn * (Md / Mdn); % Scale the data back to the original MAD
Xn = Xn + M; % Move the data to the original Median
if any(isnan(Xn))
    warning('Normalization Failed! The values were not normalized.')
    Xn = X;
    return
end


%% resume ties
for i = 1:length(tmp_Iu)
    F = tmpI==i;
    Xn(F) = median(Xn(F));
end
end