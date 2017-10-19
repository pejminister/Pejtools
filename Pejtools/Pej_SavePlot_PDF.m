function Pej_SavePlot_PDF(Fig, OutFig, Transparent)
if nargin <3; Transparent = true; end

[pathstr,~,~] = fileparts(OutFig);
if ~isempty(pathstr) && ~exist(pathstr, 'dir'); mkdir(pathstr);end

saveas(Fig, [OutFig '.fig'])
set(Fig, 'inverthardcopy', 'off');

try
    set(Fig, 'PaperPositionMode', 'auto');
catch err
end
if Transparent
    %     set(gca, 'color', 'w')
    set(gcf, 'color', 'none')
    
end 

set(Fig,'Units','Inches');
pos = get(Fig,'Position');
set(Fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

print(Fig,  [OutFig '.pdf'],'-dpdf','-r0');

close(Fig)
disp(['Figure saved: ' OutFig] )
end