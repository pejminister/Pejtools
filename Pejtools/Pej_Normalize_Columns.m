function X = Pej_Normalize_Columns(X,Offset )
ref = log(median(X,2));
XL = log(X);     clear X
if nargin<2
    Xd = XL - repmat(ref,1, size(XL,2));
    Filt = ~any(isnan(Xd),2);
    Offset = Pej_Median_withNaNs(Xd(Filt,:));
end
X = exp(XL - repmat(Offset, size(XL,1),1));
end