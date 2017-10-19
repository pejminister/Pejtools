% This function adds a f(x)=x line to the plot
function Pej_Add_IdentityDiag(Color)
if nargin==0
    Color = 'k';
end
hold on

x = min([xlim ylim]);
X = max([xlim ylim]);
plot([x X], [x X], 'color', Color);
xlim([x X])
ylim([x X])

end