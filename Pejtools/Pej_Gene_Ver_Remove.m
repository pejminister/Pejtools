% this function removes everything after a dot in the gene name:
% ENSG00000213073.4  -> ENSG00000213073
% Pej 2017, NYGC
function S = Pej_Gene_Ver_Remove(S)

for i = length(S(:)):-1:1
    t = S{i};
    f = find(t=='.', 1, 'first');
    if ~isempty(f)
        S{i}= S{i}(1:f-1);
    end
end
end