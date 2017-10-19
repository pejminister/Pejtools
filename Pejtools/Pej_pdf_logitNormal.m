% this is PDF for logit-normal distribution
% Pejman

function px = Pej_pdf_logitNormal(x,mu,Std)
px = nan(size(x));
F = x<=1 & x>=0;
px(~F) = 0;

xF = x(F);
C1 = 2*Std.^2;
px(F) = (1./(xF.*(1-xF).*sqrt(C1*pi))).*exp(-((log(xF./(1-xF))-mu).^2)./C1);

end