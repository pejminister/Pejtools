% This gets two vectors and collpses the entries in first one to a list of
% arrays associated with the unique entries in the second input (KeyField)

function [List, uniqueIDs] = Pej_Collapse(X, KeyField)
[uniqueIDs,~,I] = unique(KeyField);

for i = length(uniqueIDs):-1:1
    List(i,1)={X(I==i)};
end
end