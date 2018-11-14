% This  file joins two structures that are made by "Pej_Read_Expression_Table.m"
% either by intersecting rows and columns or by union
% Pejman, July 2018
% -------

function Joint = Pej_Expression_Table_Join(X, Y, ByUnion)
if nargin<3
    ByUnion = false;
end

if ByUnion == true
    jID = union(X.GeneNames, Y.GeneNames);
    jLbl = union(X.SampleLabels, Y.SampleLabels);
else
    jID = intersect(X.GeneNames, Y.GeneNames);
    jLbl = intersect(X.SampleLabels, Y.SampleLabels);
end

%% Match rows
[~, jia, ai] = intersect(jID, X.GeneNames);
[~, jib, bi] = intersect(jID, Y.GeneNames);

X_tmp = nan(length(jID), length(X.SampleLabels));
Y_tmp = nan(length(jID), length(Y.SampleLabels));

X_tmp(jia,:) = X.Expressions(ai,:);
Y_tmp(jib,:) = Y.Expressions(bi,:);

%% Match columns
[~, jja, aj] = intersect(jLbl, X.SampleLabels);
[~, jjb, bj] = intersect(jLbl, Y.SampleLabels);

Joint.GeneNames = jID;

Joint.Expressions1 = nan(length(jID), length(jLbl));
Joint.Expressions2 = nan(length(jID), length(jLbl));

Joint.Expressions1(:,jja) = X_tmp(:,aj);
Joint.Expressions2(:,jjb) = Y_tmp(:,bj);

Joint.SampleLabels = jLbl;
end