function S = Pej_MatlabTable2PejTable(T)
F = T.Properties.VariableNames;

for f= 1:length(F)
   S.(F{f})=T.(F{f}); 
end
end