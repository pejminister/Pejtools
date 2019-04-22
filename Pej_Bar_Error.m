% CLH(:,1) is the middle point
% CLH(:,2) is the lower  point
% CLH(:,3) is the upper  point

function Pej_Bar_Error(CLH, Labels)
if nargin < 2
    Labels = 1:size(CLH,1);
elseif length(Labels)~=size(CLH,1)
    error('Labels do not match the data')
end

bar(CLH(:,1), 'k');
hold on
errorbar(1:length(Labels), CLH(:,1), CLH(:,1)-CLH(:,2)  , nan(size(CLH(:,1))), '.', 'linewidth', 1.5, 'color', 'w', 'Marker', 'none')
errorbar(1:length(Labels), CLH(:,1), nan(size(CLH(:,1))), CLH(:,3)-CLH(:,1)  , '.', 'linewidth', 1.5, 'color', 'k', 'Marker', 'none')
set(gca, 'XTick', 1:length(Labels), 'XTicklabel', Labels, 'XTickLabelRotation',90)
xlim([.4 length(Labels)+.6])
set(gcf, 'position', [1 1 205  50*length(Labels)]);
end