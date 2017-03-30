function  Fig = Pej_Heatmap(Data, RowNames, ColumnNames, ColorRagne, ColorBarLabel)
NaNcolor = [1 1 1]*.5; % missing points
N = size(Data,1);
P = size(Data,2);
if nargin<2 || isempty(RowNames)
    RowNames = [];
end

if nargin<3 || isempty(ColumnNames)
    ColumnNames = [];
end

Fig = figure;

h = imagesc(Data);
set(h,'alphadata',~isnan(Data));
set(gca,'Color',NaNcolor);
if ~isempty(RowNames)
    if length(RowNames)~=N
        error('Row labels are different size from data')
    end
    set(gca, 'YTick', 1:N);
    set(gca, 'YTickLabel', RowNames);
    set(gca, 'YDir', 'normal');
end


if ~isempty(ColumnNames)
    if length(ColumnNames)~=P
        error('Column labels are different size from data')
    end
    set(gca, 'XTick', 1:P);
    set(gca, 'XTickLabel', ColumnNames);
try
    set(gca, 'XTickLabelRotation', 90)
catch
    
end
end

if nargin>=4 && ~isempty(ColorRagne)
   caxis(ColorRagne);
end

h = colorbar;
if nargin>=5 && ~isempty(ColorBarLabel)
    ylabel(h, ColorBarLabel);
end
end