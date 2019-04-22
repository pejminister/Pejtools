% this is for reformatting bar plots. Run it right after plotting
function Pej_Polish_Bar(h, FaceColor)
if nargin==0 || isempty(h)
    h = findobj(gca,'Type','patch');
end

if nargin<2 || isempty(FaceColor)
    FaceColor=[0 0.4470 0.7410]*.5;
end

set(h, 'FaceColor',FaceColor);
set(h, 'EdgeColor',[1 1 1]*.6);
set(h, 'linewidth',.5);
end