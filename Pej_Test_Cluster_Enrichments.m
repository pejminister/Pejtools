function Pej_Test_Cluster_Enrichments(ClusteringFile,DB_Path)
Shuffle = false; % If you put this on true, it shuffles the DEG qvalues, so should technically give flat pvalues all the time.


if nargin < 2
    DB_Path = '/Users/pmohammadi/Desktop/LocalTMP/PEJ_Resources/GeneSets/';
end
disp(['Analysing ' ClusteringFile])
Test_DEGs_Enrichments_of(ClusteringFile,DB_Path, Shuffle);

end


function Test_DEGs_Enrichments_of(Fname,DB_Path, Shuffle)
[Fldr, Fname, ~] =fileparts(Fname);

OutFldr = [ Fldr '/' Fname '_EnrichmentTests'];
mkdir(OutFldr)
Clustering_Output = [OutFldr '/Clusters_' Fname '.mat'];
tmpCluster = Pej_Read_Expression_Table([Fldr '/' Fname '.txt']);


%% Clustering
Clustered_Genes.IDs = strrep(tmpCluster.GeneNames, '"', ''); % the background for the enrichments
if Shuffle
    disp('Shuffle feature is not yet done!') % Shuffle
else
   Clustered_Genes.Cluster_IDs     = tmpCluster.Expressions==1;
   Clustered_Genes.EffectSize      = nan(size(Clustered_Genes.IDs));
   Clustered_Genes.ClusterLabels2 = tmpCluster.SampleLabels;
end


save(Clustering_Output, 'Clustered_Genes');

Pej_Test4Enrichment(Clustering_Output, DB_Path);
clear Clustered_Genes IDXKM
end