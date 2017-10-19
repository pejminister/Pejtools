% This scripts saves the figure, Fig, in the address, OutFig, in .eps, and
% .fig formats and closes it.

function Pej_SavePlot_Image(Fig, OutFig, Transparent)
if nargin <3; Transparent = true; end

[pathstr,~,~] = fileparts(OutFig);
if ~isempty(pathstr) && ~exist(pathstr, 'dir'); mkdir(pathstr);end

saveas(Fig, [OutFig '.fig'])

try
    set(Fig, 'inverthardcopy', 'off');
catch err
end
try
    set(Fig, 'PaperPositionMode', 'auto');
catch err
end
if Transparent
        set(gca, 'color', 'w')
%     set(gcf, 'color', 'none')  
end
print(Fig, [OutFig], '-dpng', '-r0', '-loose');

% set(Fig, 'PaperPositionMode', 'manual');
% set(Fig, 'PaperUnits','centimeters');
% set(Fig, 'Units','centimeters');
% pos=get(Fig,'Position');
% set(Fig, 'PaperSize', [pos(3) pos(4)]);
% set(Fig, 'PaperPosition',[0 0 pos(3) pos(4)]);
% 
% print(Fig, '-dpdf', [OutFig '.pdf']);

close(Fig)
disp(['Figure saved: ' OutFig] )
end