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
outFile = [OutFolder '/Explore_DEG_Report_' GeneSet '.txt'];
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

nFilt = isnan(Joint{1}.padj) & isnan(Joint{2}.padj); % those that were not tested in either of the DEG files
Joint{1} = Pej_Struct_RowDel(Joint{1}, nFilt);
Joint{2} = Pej_Struct_RowDel(Joint{2}, nFilt);

pFilt = ~(Joint{1}.padj<=Qthr | Joint{2}.padj<=Qthr);
cFilt = abs(Joint{1}.log2FoldChange)<Cthr & abs(Joint{2}.log2FoldChange)<Cthr;

% % Exclude those that are not significant or have much fold change
% Joint{1} = Pej_Struct_RowDel(Joint{1}, pFilt | cFilt);
% Joint{2} = Pej_Struct_RowDel(Joint{2}, pFilt | cFilt);

Difference = Joint{1}.log2FoldChange-Joint{2}.log2FoldChange;
Discodance = abs(Difference) ./ sqrt(Joint{1}.lfcSE.^2+Joint{2}.lfcSE.^2);
Pdiscord   = normcdf(Discodance, 'upper')*2;
Pdiscord(pFilt | cFilt) = nan; % Exclude those that are not significant or have much fold change

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
Report.DEG_p_1  = log10(Joint{1}.pvalue);
Report.DEG_p_2  = log10(Joint{2}.pvalue);
Report.DEG_q_1  = log10(Joint{1}.padj);
Report.DEG_q_2  = log10(Joint{2}.padj);

[~, I] = sort(Report.Discordance_p, 'ascend');
Report = Pej_Struct_RowSelect(Report, I);
Report.Index(:,1)= 1:length(Report.DEG_q_2);

Pej_Write_Table(outFile, Report);
%% Make a clustering report file, and do enrichment testing
% tqx = Joint{1}.padj;
% tqy = Joint{2}.padj;
tx  = Report.log2FoldChange_1; tx(Report.DEG_q_1>log10(Qthr))=0;
ty  = Report.log2FoldChange_2; ty(Report.DEG_q_2>log10(Qthr))=0;
tt = ~isnan(Report.Discordance_q);
td = Report.Discordance_q <= log10(Qthr);

% zones
tz(:,1) = td & tx> 0 & ty> 0 & tx< ty;
tz(:,2) = td & tx> 0 & ty> 0 & tx> ty;
tz(:,3) = td & tx> 0 & ty==0;
tz(:,4) = td & tx> 0 & ty< 0;

tz(:,5) = td & tx==0 & ty> 0;
tz(:,6) = td & tx==0 & ty==0;
tz(:,7) = td & tx==0 & ty< 0;

tz(:,8) = td & tx< 0 & ty> 0;
tz(:,9) = td & tx< 0 & ty==0;
tz(:,10)= td & tx< 0 & ty< 0 & tx> ty;
tz(:,11)= td & tx< 0 & ty< 0 & tx< ty;

tz(:,12)= ~td & tt & tx> 0;
tz(:,13)= ~td & tt & tx< 0;

Zoning = -1*ones(size(tx));
for i = 1:size(tz, 2)
    Zoning(tz(:,i)) = i;
end

Report.Zoning = Zoning;
Pej_Write_Table([outFile(1:end-4) '_zoned.txt'], Report);

% % Plot zones
% figure;
% plx = [zeros(size(tz, 2)+1,1);Joint{1}.log2FoldChange];plx(Zoning==-1)=nan;
% ply = [zeros(size(tz, 2)+1,1);Joint{2}.log2FoldChange];
% plz = [-1; [1:size(tz,2)]';Zoning];
% Pej_Scatter_Groups(plx, ply, plz);
% figure
% Pej_Hist_Percent(Zoning, -1:size(tz, 2));

figure
set(gcf, 'position', [440        -157        1687         955])
hold on
LIM = 5;

plx = Pej_Bound_to(Report.log2FoldChange_1, -LIM, LIM);
ply = Pej_Bound_to(Report.log2FoldChange_2, -LIM, LIM);

TT = [-1 1:13];
for i = 1:14
    subplot(3,5,i); hold on
    
    plot([0 0], [-1 1]* LIM, 'k')
    plot([-1 1]* LIM, [0 0], 'k')
    F = TT(i)==Zoning;
    plot(plx(F), ply(F), '.')
    title(num2str(TT(i)));
    Pej_Add_IdentityDiag;
    box on
    xlim([-1 1]*5)
    ylim([-1 1]*5)
end

subplot(3,5,15);hold on
plot([0 0], [-1 1]* LIM, 'k')
plot([-1 1]* LIM, [0 0], 'k')
F = Zoning~=-1;
Pej_Scatter_Groups(plx(F), ply(F), Zoning(F));
legend off
Pej_Add_IdentityDiag
title('All zones')
Pej_SavePlot_Image(gcf, [OutFolder '/Figures/' GeneSet '_Zoning'])


%% Enrichment test
% WARNING: THE ORDER OF THE NEXT LINES ARE IMPORTANT< DONT SHUFFLE THEM!
tmpCl.GeneNames = Report.GeneName;
tmpCl.Discordant = td;     
tmpCl.Concordant = tz(:,12) | tz(:,13);
tmpCl.U  = tz(:,12);
tmpCl.D  = tz(:,13);
tmpCl.uU = tz(:,1) | tz(:,5);
tmpCl.Uu = tz(:,2) | tz(:,3);
tmpCl.UD = tz(:,4);
tmpCl.dD = tz(:,7) | tz(:,10);
tmpCl.DU = tz(:,8);
tmpCl.Dd = tz(:,9) | tz(:,11);
Pej_Write_Table([outFile(1:end-4) '_classes.cr'], tmpCl);
Pej_Test_Cluster_Enrichments([outFile(1:end-4) '_classes.cr']);

%% plot Classes
ClassLs = fieldnames(tmpCl);
Report.ClassName = repmat({'Other'}, size(Report.Zoning,1),1);

figure
set(gcf, 'position', [440        -157        1687         600])
hold on
LIM = 5;

plx = Pej_Bound_to(Report.log2FoldChange_1, -LIM, LIM);
ply = Pej_Bound_to(Report.log2FoldChange_2, -LIM, LIM);

TT = [1:10];

for i = 1:10
    subplot(2,5,i); hold on
    plot([0 0], [-1 1]* LIM, 'k')
    plot([-1 1]* LIM, [0 0], 'k')
    F = tmpCl.(ClassLs{i+1});
    if i>3 && sum(F)>1
        % WARNIGN: shitty programming! Here I assume the order of
        % classlabels is the way variable "tmpCl" is made above, and the
        % first 3 fields are genenames, Discordant, and concordant, which
        % I'm not interetsed in. If you are changing this code you need to
        % make sure this part is not wrong!
        Report.ClassName(F) = repmat(ClassLs(i+1), sum(F),1);
    end
    plot(plx(F), ply(F), '.')
    title([ClassLs{i+1} ' (n=' int2str(sum(F)) ')']);
    Pej_Add_IdentityDiag;
    box on
    xlim([-1 1]*5)
    ylim([-1 1]*5)
end

Pej_SavePlot_Image(gcf, [OutFolder '/Figures/' GeneSet '_Classes'])

Pej_Write_Table([outFile(1:end-4) '_zoned.txt'], Report);


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
set(gcf, 'position', [  707   153   337   274]);
box on
Pej_SavePlot(Fig, [OutFolder '/Figures/' GeneSet]);

%% Add index numbers on genes
Fig = open( [OutFolder '/Figures/' GeneSet '.fig']);
set(gcf, 'position', [440   371   560   420]);
Sigs = Report.Discordance_q <= log10(Qthr) & Report.Difference_in_lfc > 0;
text(Report.log2FoldChange_1(Sigs), Report.log2FoldChange_2(Sigs), num2str(Report.Index(Sigs)), 'color', 'r', 'fontsize', 5, 'horizontalalignment', 'center', 'verticalalignment', 'middle')

Sigs = Report.Discordance_q <= log10(Qthr) & Report.Difference_in_lfc < 0;
text(Report.log2FoldChange_1(Sigs), Report.log2FoldChange_2(Sigs), num2str(Report.Index(Sigs)), 'color', 'r', 'fontsize', 5, 'horizontalalignment', 'center', 'verticalalignment', 'middle')

Pej_SavePlot(Fig, [OutFolder '/Figures/' GeneSet '_Indexed']);

end

function DEG_Res = Read_DEG(DEGoutputfile)
DEG_Res = Pej_Read_Table(DEGoutputfile);
% DEG_Res = Pej_Struct_RowDel(DEG_Res,isnan(DEG_Res.padj));
DEG_Res.GeneIDs = strrep(DEG_Res.UnLabeled_C1, '"', '');

tmpID = Pej_Xref(DEG_Res.GeneIDs, '/Users/pmohammadi/Desktop/LocalTMP/PEJ_Resources/Gene-Xref-25Apr2015.txt');
DEG_Res.GeneNames = tmpID.Associated_Gene_Name;
end