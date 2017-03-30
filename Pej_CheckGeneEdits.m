% This gets one folder, and in there looks for fastq files. It loads the
% WHOLE fatsq in RAM and compares it to a reference (or its reverse complement)
% which is provded in the same folder with the same name as the fastqfile
% like:
% myfile.fastq
% myfile_ref.fasta
% then for each file it reports the observed substitutions and indels up to
% a point.
% this program only looks into the area of the reference locus that at least 50 percent of the reads are aligned to it.
% Pej May 2015


% The second way to use the code is to give it a directory and a given
% reference to be used for all the fastq files in the directory, in this
% case the previous format will be ignored. (Pej Sept 2015)

% Example:
% DataDir = '/Users/pmohammadi/Desktop/LocalTMP/Ana_GeneEdits_Apr2016/PEX6_M_April2016_Trimmed/PEX6_M_data';
% RefName = 'PEX6_M_PCR.ref.txt'; 
% TemplateName = 'PEX6_M.template';
% BatchMode = false;
% Pej_CheckGeneEdits(DataDir, RefName, TemplateName, BatchMode)

%------------


function Pej_CheckGeneEdits(DataDir, RefName, TemplateName, BatchMode)
try
    addpath('/data/research/tuuli_lab/pmohammadi/PejLib')
catch er
end

if nargin<1; DataDir='.'; end
if nargin<3; TemplateName=[]; end
if nargin<4; BatchMode=false; end

AllFiles = dir(DataDir);
if isempty(AllFiles)
    fprintf('No Files found under the provided path!! Check if the path is correct\n')
end
% if BatchMode
    SummaryFoutName= [DataDir '/Summary_Report.txt'];
    SummaryFout = fopen(SummaryFoutName, 'w' );
    fprintf(SummaryFout, ...
        'File\tFileSize(MB)\tAnalysisStatus\tReverseComplemented\tTotalReads\tDiscardedReads\tDiscardedReads%%\tEditionLocation\tReference\tTemplate\tNoIndelReads\tEditedReads\tEditedReads%%\tA%%\tC%%\tG%%\tT%%\n');
    fclose(SummaryFout);
% end

% Pej_GetFiles([DEGoutputfolder '/*.txt']);
% Resultfiles = regexp(Resultfiles, '[\f\n\r]', 'split');
for f = 1 : length(AllFiles)
    if ~isempty(AllFiles)
        [~, Fname, Fext] = fileparts(AllFiles(f).name);
        if ~strcmpi(Fext, '.fastq')
            continue
        end
        disp(['Analysing ' AllFiles(f).name])
        if nargin<2 || isempty(RefName)
            Test_Align_Final_sub([DataDir '/' Fname Fext], [DataDir '/' Fname '_ref.fasta'], [], SummaryFoutName);
        else
            Test_Align_Final_sub([DataDir '/' Fname Fext], [DataDir '/' RefName], [DataDir '/' TemplateName], SummaryFoutName);
        end
    end
end
end
function Test_Align_Final_sub(InFile, InFileRef, TemplateName, SummaryFoutName)
Reported = false;
APthr = 1E-3; % pvalue threshold of the alignment score
TotReads = Inf; % Total reads to analyze
SampleReads = 1000; % Max Number of reads sampled for evaluating Score disctributions
try
    Reads2 = fastqread(InFile, 'Blockread', [1 TotReads]);
    Locus = Get_Seq(InFileRef);
    if nargin>2 && ~isempty(TemplateName)
        
        Template = Get_Seq(TemplateName);
        
        [Ts_vs_T    , Ref_vs_Template] = nwalign(Locus,Template, 'Glocal', true);
        [Ts_vs_T_rc , Ref_vs_Template_rc] = nwalign(Locus,seqrcomplement(Template), 'Glocal', true);
        
        if Ts_vs_T_rc > Ts_vs_T
            % Reverse complement the template
            Ref_vs_Template = Ref_vs_Template_rc;
            Template = seqrcomplement(Template);
        end
        IsNotDel_in_Template = Ref_vs_Template(3,:)~='-';
        
        TIf = find(IsNotDel_in_Template,1,'first');
        TIl = find(IsNotDel_in_Template,1,'last');
        if ~all(IsNotDel_in_Template(TIf:TIl))
            error('There are deletions in the template! Not supported!')
        end
        if any(Ref_vs_Template(1,:)=='-')
            error('There are insertions in the template! Not supported!')
        end
        
        I_Mismatch_in_Template = TIf + find(Ref_vs_Template(3,TIf:TIl)~=Ref_vs_Template(1,TIf:TIl)) - 1;
        if isempty(I_Mismatch_in_Template)
            error('There are no mismatches in the template! Not supported!')
        end
        if length(I_Mismatch_in_Template)>1
            error('There are more than one mismatches in the template! Not supported!')
        end
        
        
        
        Template_Target = Ref_vs_Template(3, I_Mismatch_in_Template);
        Template_TargetRef = Ref_vs_Template(1, I_Mismatch_in_Template);
    end
    L = length(Locus);
    N = length(Reads2);
    SampleReads = min(N,SampleReads);
    Reads2 = {Reads2(:).Sequence};
    
    %% Find which way to align it!
    fprintf('Calculating the expected alignment scores...')
    t = randperm(N,SampleReads);
    tmpS = nan(SampleReads, 4);
    for i = 1:SampleReads
        LR =length(Reads2{t(i)});
        RCL = seqrcomplement(Locus);
        
        tmpS(i,1) = nwalign(Locus,Reads2(t(i)), 'Glocal', true);
        tmpS(i,3) = nwalign(RCL  ,Reads2(t(i)), 'Glocal', true);
        tmpS(i,2) = nwalign(Locus,Reads2{t(i)}(randi(LR,LR,1)), 'Glocal', true);
        tmpS(i,4) = nwalign(RCL  ,Reads2{t(i)}(randi(LR,LR,1)), 'Glocal', true);
        
    end
    fprintf('Done!\n')
    
    
    if median(tmpS(:,1))<median(tmpS(:,3))
        % reverse complement all the reads!
        for i = 1:N
            Reads2{i} = seqrcomplement(Reads2{i});
        end
        RVflag=true;
        AscoreTHR = std(tmpS(:,[4])*norminv(1-APthr)) + median(tmpS(:,[4]));
    else
        RVflag = false;
        AscoreTHR = std(tmpS(:,[2])*norminv(1-APthr)) + median(tmpS(:,[2]));
    end
    
    
    %% the rest
    tTXT=0;
    Pat = zeros(16,L+1);
    Ins_Cnt_Len=zeros(L+1,2);
    Del_Cnt_Len=zeros(L+1,2);
    Discarded=0;
    hist_nMismatch_in_HR = zeros(L,1);
    nHR = 0; % number of homologous recombinant reads
    nHR_Background = zeros(4,1); % number of reads that qualify as homologous recombinant reads except that they have a wrong nt
    Score(N)=0;Alignment_{N}='';
    parfor i = 1:N
        [Score(i), Alignment_{i}] = nwalign(Locus,Reads2{i}, 'Glocal', true);
    end
    
    for i = N:-1:1
        Alignment = Alignment_{i};
        if Score(i)>AscoreTHR
            
        else
            % discard it!
            Discarded = Discarded+1;
            continue
        end
        
        %% Check for Homologous Recombinantions (HR)
        tmpA = Alignment;
        
        [MiddleGap, SideGapFilt]= Has_Gap_in_Middle(tmpA(1,:), L);
        if MiddleGap
            % its not a HR
        else
            tmpA = tmpA(:,SideGapFilt);
            [MiddleGap, SideGapFilt2]= Has_Gap_in_Middle(tmpA(3,:), length(Reads2{i}));
            if MiddleGap
                % its not a HR
            else
                if upper(tmpA(1,I_Mismatch_in_Template)) ~= Template_TargetRef
                    % this means this code is buggy!!
                    error('ups!!!!!!')
                end
                
                if upper(tmpA(3,I_Mismatch_in_Template)) == Template_Target
                    % Bingo! Found an HR!
                    nHR = nHR + 1;
                    nmHR = sum(SideGapFilt2 & (tmpA(1,:)~=tmpA(3,:)));
                    hist_nMismatch_in_HR(nmHR) = hist_nMismatch_in_HR(nmHR)+1;
                else
                    switch  upper(tmpA(3,I_Mismatch_in_Template))
                        case 'A'
                            nHR_Background(1) = nHR_Background(1) + 1;
                        case 'C'
                            nHR_Background(2) = nHR_Background(2) + 1;
                        case 'G'
                            nHR_Background(3) = nHR_Background(3) + 1;
                        case 'T'
                            nHR_Background(4) = nHR_Background(4) + 1;
                    end
                    
                end
                
            end
        end
        
        
        
        %% Get patterns
        tmpA = Alignment;
        tmpA(:,tmpA(1,:)=='-')=[];
        
        tmpint = nt2int(tmpA(3,:));
        C = unique(tmpint);
        
        g=1;
        while tmpA(3, g)=='-' && g<L
            % remove deletions at the begining
            tmpint(g)=nan;
            g=g+1;
        end
        g=L;
        while tmpA(3, g)=='-' && g>1
            % remove deletions at the end
            tmpint(g)=nan;
            g=g-1;
        end
        
        for j =1:length(C)
            Fi=tmpint==C(j);
            Pat(C(j),Fi)= Pat(C(j),Fi)+1;
        end
        
        
        
        tmpA = Alignment([1 3],:);
        %% Number of reads with a insertions
        [tGaps_Starts, tGaps_Length] = Get_Gaps(tmpA(1,:));
        for j = 1:length(tGaps_Starts)
            tPos = tGaps_Starts(j) - sum(tGaps_Length(1:j-1));
            Ins_Cnt_Len(tPos,1)= Ins_Cnt_Len(tPos,1)+1;
            Ins_Cnt_Len(tPos,2)= Ins_Cnt_Len(tPos,2)+tGaps_Length(j);
        end
        
        %% Number of reads with a deletions
        tmpA(:,tmpA(1,:)=='-')=[];
        [tGaps_Starts, tGaps_Length] = Get_Gaps(tmpA(2,:));
        if ~isempty(tGaps_Starts) && tGaps_Starts(1)==1
            % it the begining of the string discard it
            tGaps_Starts(1)=[];tGaps_Length(1)=[];
        end
        
        if ~isempty(tGaps_Starts) && tGaps_Starts(end)+tGaps_Length(end)==L+1
            % it the end of the string discard it
            tGaps_Starts(end)=[];tGaps_Length(end)=[];
        end
        
        Del_Cnt_Len(tGaps_Starts,1)= Del_Cnt_Len(tGaps_Starts,1)+1;
        Del_Cnt_Len(tGaps_Starts,2)= Del_Cnt_Len(tGaps_Starts,2)+tGaps_Length';
        
        if (i/1000)==round(i/1000)
            fprintf(repmat('\b', 1, tTXT));
            tTXT = fprintf('%d out of %d', i, N);
        end
        
    end
    
    %% Do all detailed report stuff
    
    figure
    % boxplot(tmpS);
    Pej_Plot_Dist_Custom([1:4], tmpS, [], 1)
    Pej_Plot_Dist_Custom([1:4], tmpS, [], -1)
    
    set(gca, 'XTick', 1:4)
    set(gca, 'xticklabel', {'-', '-*', 'RC', 'RC*'})
    ylabel('Alignment Score')
    box on
    xlim([.5 4.5])
    set(gcf, 'position', [1 1 400 300]);
    plot([.8 4.2], AscoreTHR*[1 1], '--', 'color', [.8 .1 .05]);
    
    if RVflag
        text(3, AscoreTHR, sprintf('%d(%.2f%%)\ndiscarded', Discarded, Discarded/N*100), 'horizontalalignment', 'left', 'verticalalignment', 'top', 'color', [.6 0 0]);
        errorbar([ 4], median(tmpS(:,[ 4])),[0 ], std(tmpS(:,[4])*norminv(1-APthr)), '.', 'linewidth', 1.5, 'color', [.8 .1 .05]);
    else
        text(1, AscoreTHR, sprintf('%d(%.2f%%)\ndiscarded', Discarded, Discarded/N*100), 'horizontalalignment', 'left', 'verticalalignment', 'top', 'color', [.6 0 0]);
        errorbar([2], median(tmpS(:,[2])),[0], std(tmpS(:,[2 ])*norminv(1-APthr)), '.', 'linewidth', 1.5, 'color', [.8 .1 .05]);
        
    end
    Pej_SavePlot(gcf,[InFile '_results/AlignmentQC']);
    
    
    mX = find(sum(Pat,1)>.10*N, 1, 'first');
    if isempty(mX)
        mX = 1;
    end
    MX = find(sum(Pat,1)>.10*N, 1, 'last');
    MX = max(MX, mX+10);
    
    
    F=Del_Cnt_Len(:,1)>0;
    Del_Cnt_Len(F,2) = Del_Cnt_Len(F,2)./Del_Cnt_Len(F,1);
    
    F=Ins_Cnt_Len(:,1)>0;
    Ins_Cnt_Len(F,2) = Ins_Cnt_Len(F,2)./Ins_Cnt_Len(F,1);
    
    if RVflag
        Fout = fopen([InFile '_results/Pileup_' int2str(N) '_ReadsRCed.txt'], 'w');
    else
        Fout = fopen([InFile '_results/Pileup_' int2str(N) '_Reads.txt'], 'w');
    end
    
    Fpat = sum(Pat,2)==0;
    Patguide = 1:16;
    
    Pat(Fpat,:)=[];
    Patguide(Fpat)=[];
    TotCov = sum(Pat,1);
    
    outBuff=[[1:L+1]; double([Locus '-']); TotCov; Del_Cnt_Len'; Ins_Cnt_Len'; Pat];
    
    fprintf(Fout, '# Template is different from Reference in position %d (%c=>%c)\n',I_Mismatch_in_Template,  Template_TargetRef, Template_Target);
    fprintf(Fout, '# Percentage of Homologous Recombinations: %.2f%% (%d/%d)\n', nHR/TotCov(I_Mismatch_in_Template)*100, nHR, TotCov(I_Mismatch_in_Template));
    fprintf(Fout, '# Background percentage of Homologous Recombinations:\tA:%.2f%% C:%.2f%% G:%.2f%% T:%.2f%% \n', nHR_Background/TotCov(I_Mismatch_in_Template)*100);
    
    fprintf(Fout, ['Position\tReference\tTotalCoverage\tDelCount\tDelLength\tInsCount\tInsLength' repmat('\t%c', 1, size(Pat,1)) '\n'],int2nt(Patguide));
    fprintf(Fout, ['%d\t%c\t%d\t%d\t%.1f\t%d\t%.1f', repmat('\t%d', 1, size(Pat,1)) '\n'],outBuff);
    fclose(Fout);
    
    %% Make Batch report
    SummaryFout = fopen(SummaryFoutName, 'a' );
    InFileInfo = dir(InFile);
    InFileSizeMB = InFileInfo.bytes/1E+6;
    
    %'File\tFileSize(MB)\tAnalysisStatus\tReverseComplemented\tTotalReads\tDiscardedReads\tDiscardedReads%%\tEditionLocation\tReference\tTemplate\tNoIndelReads\tEditedReads\tEditedReads%%\tA%%\tC%%\tG%%\tT%%\n');
    fprintf(SummaryFout, '%s\t%.3f\t%s\t%d\t%d\t%d\t%.2f\t%d\t%c\t%c',InFile,InFileSizeMB, 'Worked', RVflag, N, Discarded, Discarded/N*100,I_Mismatch_in_Template,  Template_TargetRef, Template_Target);
    fprintf(SummaryFout, '\t%d\t%d\t%.3f', TotCov(I_Mismatch_in_Template), nHR, nHR/TotCov(I_Mismatch_in_Template)*100);
    fprintf(SummaryFout, '\t%.3f\t%.3f\t%.3f\t%.3f\n', nHR_Background/TotCov(I_Mismatch_in_Template)*100);
    fclose(SummaryFout);
    Reported = true;
    %% HR plot
    figure;
    tM = find(hist_nMismatch_in_HR>0,1, 'last');
    if isempty(tM) || tM<5
        tM =5;
    end
    bar(1:tM, hist_nMismatch_in_HR(1:tM)/nHR*100, 'hist');
    
    ylabel('% of HR reads')
    xlabel('Number of mismatches in HR reads')
    Pej_SavePlot(gcf, [InFile '_results/HRecomb_Mismatch_Dist']);
    
    %% Next
    figure;
    t = Del_Cnt_Len(:,1)+Ins_Cnt_Len(:,1);
    Borders(1) = min(find(t, 1, 'first'));
    Borders(2) = min(find(t, 1, 'last'));
    Borders(1) = floor(Borders(1)/10)*10;
    Borders(2) = floor(Borders(2)/10)*10;
    
    subplot(2,1,1)
    Del_Cnt_Len(Del_Cnt_Len(:,1)==0,1)=1E-1;
    bar(1+log10(Del_Cnt_Len(:,1)));hold on
    bar(-Del_Cnt_Len(:,2));
    
    YT = get(gca, 'YTick');
    c=1;
    while YT(c)<0
        c=c+1;
        YTL{c}=sprintf('%d', abs(YT(c)));
    end
    while c<length(YT)
        c=c+1;
        YTL{c}=sprintf('10^{%.1f}', YT(c)-1);
    end
    set(gca, 'YTick', YT);
    set(gca, 'YTickLabel', YTL);
    
    grid on
    xlabel('Deletion at position on reference locus (nt)')
    ylabel('Average length / Count')
    xlim(Borders)
    
    subplot(2,1,2)
    Ins_Cnt_Len(Ins_Cnt_Len(:,1)==0,1)=1E-1;
    bar(1+log10(Ins_Cnt_Len(:,1)));hold on
    bar(-Ins_Cnt_Len(:,2));
    xlim(Borders)
    
    grid on
    xlabel('Insertion before position on reference locus (nt)')
    ylabel('Average length / Count')
    
    YT = get(gca, 'YTick');
    c=1;
    while YT(c)<0
        c=c+1;
        YTL{c}=sprintf('%d', abs(YT(c)));
    end
    while c<length(YT)
        c=c+1;
        YTL{c}=sprintf('10^{%.1f}', YT(c)-1);
    end
    set(gca, 'YTick', YT);
    set(gca, 'YTickLabel', YTL);
    
    if RVflag
        Pej_SavePlot(gcf, [InFile '_results/Indels_' int2str(N) '_ReadsRCed']);
    else
        Pej_SavePlot(gcf, [InFile '_results/Indels_' int2str(N) '_Reads']);
    end
    
    
    
    
    
    
    % %%
    
    
    Pat(:,1:mX)=nan;
    Pat(:,MX:end)=nan;
    
    
    Pat = Pat./repmat(TotCov,size(Pat,1),1)*100;
    
    
    inLocus = nt2int([Locus '-']);
    for l = 1:size(Pat,2)
        Fl = Patguide==inLocus(l);
        Pat(Fl,l)=0;
    end
    figure;
    colormap(jet);
    bar(Pat', 'stacked')
    legend(int2nt(Patguide'))
    xlabel('Position on reference locus (nt)')
    ylabel('Mismatch ratio(%)')
    set(gcf, 'position',  [1634         673        2000         200])
    set(gca, 'position', [0.05    0.2    0.9    0.75])
    ylim([0, max(ylim)+.5]);
    xlim([mX, MX])
    set(gca, 'color', 'none');
    grid on
    
    if RVflag
        Pej_SavePlot(gcf, [InFile '_results/Mismatches_' int2str(N) '_ReadsRCed']);
    else
        Pej_SavePlot(gcf, [InFile '_results/Mismatches_' int2str(N) '_Reads']);
    end
catch err
    err.getReport
    try
        SummaryFout = fopen(SummaryFoutName, 'a' );
        InFileInfo = dir(InFile);
        InFileSizeMB = InFileInfo.bytes/1E+6;
    catch
        InFileSizeMB = nan;
    end
    try
        if ~Reported
            SummaryFout = fopen(SummaryFoutName, 'a' );
            %'File\tFileSize(MB)\tAnalysisStatus\tReverseComplemented\tTotalReads\tDiscardedReads\tDiscardedReads%%\tEditionLocation\tReference\tTemplate\tNoIndelReads\tEditedReads\tEditedReads%%\tA%%\tC%%\tG%%\tT%%\n');
            fprintf(SummaryFout, '%s\t%.3f\t%s\t%d\t%d\t%d\t%c\t%c',InFile, InFileSizeMB, 'Failed!', nan, nan,nan,  '-', '-');
            fprintf(SummaryFout, '\t%d\t%d\t%.3f', nan, nan, nan);
            fprintf(SummaryFout, '\t%.3f\t%.3f\t%.3f\t%.3f\n', nan(4,1));
            fclose(SummaryFout);
        end
        close all
        fclose all
        
    catch
        warning('Really failed!!!')
    end
end

end


function Pej_Plot_Dist_Custom(x,data, Color, Direction, Dashed, Scale)
if nargin<2; data = x; x = 1:size(data,2); end

MidPoints =  median(data);

if strcmpi(get(gca, 'YScale'), 'log')
    %% estimate densities from the log(data), and plot it back on a semilogY scale.
    data = log(data);
    isLog = true;
else
    isLog = false;
end

N = size(data,2);
if nargin<3 || isempty(Color);       Color      = repmat([.05,0.1,.2], N, 1); end
if nargin<4 || isempty(Direction);   Direction  = ones(N,1); end
if nargin<5 || isempty(Dashed);      Dashed     = false; end
if nargin<6 || isempty(Scale);       Scale      = nan; end

if size(Color,1)    == 1; Color     = repmat(Color    ,N,1); end
if size(Direction,1)== 1; Direction = repmat(Direction,N,1); end
if length(Scale)    == 1; Scale     = repmat(Scale    ,N,1); end

for i = 1:N
    [f,xi] = ksdensity(data(:,i),linspace(min(data(:,i)), max(data(:,i)), 300));
    if isnan(Scale(i))
        f = f/max(f) * .4;
    else
        f = f * Scale(i);
    end
    if Dashed
        
        tmpF = [f; f];
        tmpF(1,1:2:end)=0; tmpF(2,2:2:end)=0;
        tmpXi= [xi;xi];
        
        f = tmpF(:)';
        xi = tmpXi(:)';
    end
    
    X975 = quantile(data(:,i), 1);
    X025 = quantile(data(:,i), 0);
    f(xi>X975)  = []; f(xi<X025)  = [];
    xi(xi>X975) = []; xi(xi<X025) = [];
    if isLog
        xi = exp(xi);
    end
    hold on
    
    fill(x(i) + Direction(i) * [0 f 0 0],[xi(1) xi+sqrt(eps) xi(end) xi(1)], Color(i,:), 'edgecolor', Color(i,:));
    %     plot(x(i), MidPoints(i), 'o', 'color', [1 1 1], 'linewidth', .5   , 'markersize', 6, 'markerfacecolor', 0*[1 1 1])
    %     plot(x(i), MidPoints(i), '+', 'color', [1 1 1], 'linewidth', 1  , 'markersize', 4 )
end
end

function [tGaps_Starts, tGaps_Length] = Get_Gaps(S)

tmpA_Ins    = diff(['*' S '*'] =='-');
Flips = find(tmpA_Ins);

tGaps_Starts = zeros(1,length(Flips));
tGaps_Length = zeros(1,length(Flips));

if ~isempty(tGaps_Starts)
    tGaps_Length = Flips(2:2:end) - Flips(1:2:end);
    tGaps_Starts = Flips(1:2:end);
end

end

function Locus = Get_Seq(InFileRef)
Fid = fopen(InFileRef, 'r');
Locus= fgetl(Fid);
if Locus(1)=='>'
    % it's fasta formatted, read the 2nd line
    Locus= fgetl(Fid);
end
fclose(Fid);
Locus = upper(Locus);
Locus(Locus==' ')=[];
end


function [MiddleGap, SideGapFilt]= Has_Gap_in_Middle(tmpA, ReadLength)
I1 = -1; % cut from begining till here
I2 = 1+length(tmpA); % cut from here till end
MiddleGap = false;
SideGapFilt = true(size(tmpA));

[Inserts_b, Inserts_l] = Get_Gaps(tmpA);
if ~isempty(Inserts_b)
    if length(Inserts_b)>2
        % Read is not a HR
        MiddleGap = true;
    else
        % There's an instertion
        if Inserts_b(1)==1
            % insertion is in the begining
            I1 = Inserts_l(1);
            tmpA(1:I1)=[];
            [Inserts_b, ~] = Get_Gaps(tmpA);
        end
        
        if isempty(Inserts_b)
            % There's no more insertions
        else
            if Inserts_b(1)==ReadLength+1
                % Insert is at the end of the read
                I2 = ReadLength+I1+1;
            else
                % Read is not a HR
                MiddleGap = true;
            end
        end
    end
end

SideGapFilt(1:I1)= false;
SideGapFilt(I2:end)= false;
end