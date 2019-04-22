% CellLabel is the same size as "Data", and it can be cell array or a numeric array.
function  Fig = Pej_Heatmap(Data, RowNames, ColumnNames, ColorRagne, ColorBarLabel, CellLabel)
NaNcolor = [1 1 1]*.5; % missing points
N = size(Data,1);
P = size(Data,2);
if nargin<2 || isempty(RowNames)
    RowNames = [];
end

if nargin<3 || isempty(ColumnNames)
    ColumnNames = [];
end

if nargin<6 || isempty(CellLabel)
    CellLabel = [];
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
    set(gca, 'YDir', 'reverse');
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

if nargin>=5 && ~isempty(CellLabel)
    
    for i = 1: N
        for j = 1: P
            if isnumeric(CellLabel)
                if ~isnan(CellLabel(i,j))
                    text(j , i, num2str(CellLabel(i,j)), 'verticalalignment', 'middle', 'horizontalalignment', 'center');
                end
            else
                text(j , i, CellLabel{i,j}, 'verticalalignment', 'middle', 'horizontalalignment', 'center');
            end
        end
    end
end

s(1) = P/N; s(2)=1;
s= s/(min(s));

set(gcf, 'position', [100 100 200*s(1)+20 200*s(2)])
end