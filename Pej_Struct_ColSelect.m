% This selects the columns that are names in "ColNames" out of "S"
% example: 
function Sc = Pej_Struct_ColSelect(S, ColNames)
Sc = rmfield(S, setdiff(fieldnames(S), ColNames));
end