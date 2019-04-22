% This function tells you the absolute size of the active plot within the active figure

 function [W2Hratio, Width, Height] = Pej_Plot_AbsSize()

t = get(gcf,'position').*get(gca, 'position');
Width = t(3);
Height= t(4);

W2Hratio = Width/Height;
end