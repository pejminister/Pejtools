% This gets astructure and collpses the entries over one of the columns,
% using union. For example you can use this to collape info from
% a bunch of SNPs based on Gene Name. 

function Sc = Pej_Struct_Collapse(S, KeyField)
[~,I] = sort(S.(KeyField));
S = Pej_Struct_RowSelect(S,I); % sort by the key field

uL = length(unique(S.(KeyField)));

FN = fieldnames(S);
Curr = ''; i0=length(I); iu=uL;

for i = length(I):-1:1
    if isequal(Curr, S.(KeyField)(i)) %&& ~strcmp('NaN', Curr)
        % same gene
    
    else
        % New gene, summarize the previous gene and write it
        Ui=i0:i+1;
        for j = 1:length(FN)
            Sc.(FN{j}) = union()
        end 
end
end