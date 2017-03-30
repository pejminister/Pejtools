% This gets two DESeq2 result fiiles and compares them within the context
% of some Genesets. It uses the logfoldchange stds, for evaluating the
% significance of the deviation in the two DEGs.

% Example:
% Pej_Explore_DEG_GeneSet(...
%     '/Users/pmohammadi/Dropbox/HIV_MRJ_Edits/04_DEG_tests/02/DESeq_Results/DEG02.txt',...
%     '/Users/pmohammadi/Dropbox/HIV_MRJ_Edits/04_DEG_tests/03/DESeq_Results/DEG03.txt',...
%     '/Users/pmohammadi/Desktop/LocalTMP/PEJ_Resources/GeneSets/ReactomePathways-27Apr2015.gmt GeneSets.mat',...
%     {'Signal Transduction', 'Immune System'},...
%     'Stim_in_B2MvsIFI16')
%
% Pej_Explore_DEG_GeneSet(...
%     '/Users/pmohammadi/Dropbox/HIV_MRJ_Edits/04_DEG_tests/01/DESeq_Results/DEG01.txt',...
%     '/Users/pmohammadi/Dropbox/HIV_MRJ_Edits/04_DEG_tests/04/DESeq_Results/DEG04.txt',...
%     '/Users/pmohammadi/Desktop/LocalTMP/PEJ_Resources/GeneSets/ReactomePathways-27Apr2015.gmt GeneSets.mat',...
%     {'Signal Transduction', 'Immune System'},...
%     'IFI16Knockout_in_CtrlvsStim')


% Pej June 22nd, 2015,
% Chicago University, GTEx meeting

function Pej_Explore_DEG_GeneSet(DEGoutputfile1, DEGoutputfile2, ImportedGeneSetFile, GeneSetsToVisualize, OutFolder)
DEG_Res1 = Read_DEG(DEGoutputfile1);
DEG_Res2 = Read_DEG(DEGoutputfile2);

[~, F1_name, ~] = fileparts(DEGoutputfile1);
[~, F2_name, ~] = fileparts(DEGoutputfile2);
if nargin<5
    OutFolder = [F1_name '-vs-' F2_name];
    mkdir(OutFolder);
end

%% Read in DEG files
if isempty(ImportedGeneSetFile)
    % Do it for all genes, there's no geneset
    GeneSetsToVisualize = 'All_Genes';
    Pej_Explore_DEG_GeneSet_sub(DEG_Res1,DEG_Res2, [], GeneSetsToVisualize, OutFolder);
else
    %% Read in the GeneSets
    GeneSets = load(ImportedGeneSetFile);
    
    %% Run for each Geneset
    if ~iscell(GeneSetsToVisualize)
        GeneSetsToVisualize = {GeneSetsToVisualize};
    end
    for k = 1: length(GeneSetsToVisualize)
        Pej_Explore_DEG_GeneSet_sub(DEG_Res1,DEG_Res2, GeneSets, GeneSetsToVisualize{k},OutFolder);
    end
end
end

function Pej_Explore_DEG_GeneSet_sub(DEG_Res1,DEG_Res2, AllSets, GeneSet, OutFolder)
Qthr = .05;% cut off on the pvalues, both for including genes in the analysis and for the plotting the final discrepancy qvalue
Cthr = 0;%.5; % cut-off on effect size in log2 scale

if isempty(AllSets)
    % Do it for all genes, there's no geneset
else
    I =[];
    for i = 1:length(AllSets.GeneSets)
        if strcmpi(AllSets.GeneSets{i}.Name{1}, GeneSet)
            I = i;
            break
        end
    end
    if isempty(I)
        % GeneSet not found!
        warning('%s not found!', GeneSet);
        return
    end
    
    %% limit to this set
    Set = AllSets.GeneSets{I};
    Set.GeneNames = AllSets.GeneNames(Set.Genes_PejIDs);
    
    [~, ai1, bi1] = intersect(Set.GeneNames, DEG_Res1.GeneNames);
    [~, ai2, bi2] = intersect(Set.GeneNames, DEG_Res2.GeneNames);
    
    % Exclude those not in the gene set
    DEG_Res1 = Pej_Struct_RowSelect(DEG_Res1, bi1);
    DEG_Res2 = Pej_Struct_RowSelect(DEG_Res2, bi2);
end
%% Join sets
Joint = Pej_Struct_Join({'GeneIDs'}, {DEG_Res1, DEG_Res2}, true);
% Exclude those that are not significant or have much fold change
pFilt = ~(Joint{1}.padj<=Qthr | Joint{2}.padj<=Qthr);
cFilt = abs(Joint{1}.log2FoldChange)<Cthr & abs(Joint{2}.log2FoldChange)<Cthr;

Joint{1} = Pej_Struct_RowDel(Joint{1}, pFilt | cFilt);
Joint{2} = Pej_Struct_RowDel(Joint{2}, pFilt | cFilt);

Difference = Joint{1}.log2FoldChange-Joint{2}.log2FoldChange;
Discodance = abs(Difference) ./ sqrt(Joint{1}.lfcSE.^2+Joint{2}.lfcSE.^2);
Pdiscord   = normcdf(Discodance, 'upper')*2;
Qdiscord   = mafdr(Pdiscord, 'bhfdr', true);

%% Make a report file
Report.GeneIDs = Joint{1}.GeneIDs;
Report.GeneName = Joint{1}.GeneNames;
Report.Difference_in_lfc = Difference;
Report.Discordance_p = log10(Pdiscord);
Report.Discordance_q = log10(Qdiscord);
Report.log2FoldChange_1=Joint{1}.log2FoldChange;
Report.log2FoldChange_2=Joint{2}.log2FoldChange;
Report.lfcSE_1  = Joint{1}.lfcSE;
Report.lfcSE_2  = Joint{2}.lfcSE;
Report.DEG_p_1  = Joint{1}.pvalue;
Report.DEG_p_2  = Joint{2}.pvalue;
Report.DEG_q_1  = Joint{1}.padj;
Report.DEG_q_2  = Joint{2}.padj;

[~, I] = sort(Report.Discordance_p, 'ascend');
Report = Pej_Struct_RowSelect(Report, I);
Report.Index(:,1)= 1:length(Report.DEG_q_2);

Pej_Write_Table([OutFolder '/Explore_DEG_Report_' GeneSet '.txt'], Report);
%% plot it all!
Cbox = gray(256);
Cbox(end:-1:1) = Cbox;
C = -Report.Discordance_q;
MC = min(10);
C = 1+ceil(C./MC*(length(Cbox)-1));
C(C>length(Cbox))=length(Cbox);
Fig = figure;
hold on
[~, plotI] = sort(C);
for i = 1:length(Report.log2FoldChange_1)
    if ~isnan(C(plotI(i)))
        plot(Report.log2FoldChange_1(plotI(i)), Report.log2FoldChange_2(plotI(i)), 's', 'markerfacecolor', Cbox(C(plotI(i)),:), 'color',Cbox(C(plotI(i)),:)/2, 'linewidth', .1)
        if  Report.Discordance_q(plotI(i))<=log10(Qthr);
            plot(Report.log2FoldChange_1(plotI(i)), Report.log2FoldChange_2(plotI(i)), '.', 'color',[0 0 0])
        end
    end
end
colormap(Cbox);
caxis([1 MC]);
tL = colorbar;
set(get(tL,'ylabel'),'String', '-log_{10} p-value');
M = max(abs([xlim ylim]));
xlim([-M M])
ylim([-M M])
plot([0 0] , [-M M], 'k')
plot([-M M], [0  0], 'k')
plot([-M M], [-M M], '--k')
xlabel('1^{st} fold change (log_2)')
ylabel('2^{nd} fold change (log_2)')
Pej_SavePlot(Fig, [OutFolder '/Figures/' GeneSet]);

%% Add index numbers on genes
Fig = open( [OutFolder '/Figures/' GeneSet '.fig']);

Sigs = Report.Discordance_q <= log10(Qthr) & Report.Difference_in_lfc > 0;
text(Report.log2FoldChange_1(Sigs), Report.log2FoldChange_2(Sigs), num2str(Report.Index(Sigs)), 'color', 'r', 'fontsize', 5, 'horizontalalignment', 'center', 'verticalalignment', 'middle')

Sigs = Report.Discordance_q <= log10(Qthr) & Report.Difference_in_lfc < 0;
text(Report.log2FoldChange_1(Sigs), Report.log2FoldChange_2(Sigs), num2str(Report.Index(Sigs)), 'color', 'r', 'fontsize', 5, 'horizontalalignment', 'center', 'verticalalignment', 'middle')
Pej_SavePlot(Fig, [OutFolder '/Figures/' GeneSet '_Indexed']);

end

function DEG_Res = Read_DEG(DEGoutputfile)
DEG_Res = Pej_Read_Table(DEGoutputfile);
DEG_Res = Pej_Struct_RowDel(DEG_Res,isnan(DEG_Res.padj));
DEG_Res.GeneIDs = strrep(DEG_Res.UnLabeled_C1, '"', '');

tmpID = Pej_Xref(DEG_Res.GeneIDs, '/Users/pmohammadi/Desktop/LocalTMP/PEJ_Resources/Gene-Xref-25Apr2015.txt');
DEG_Res.GeneNames = tmpID.Associated_Gene_Name;
end