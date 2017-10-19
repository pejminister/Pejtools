function Pej_Write_Table(FileName, Data, ColLabels, RowLabels)
[pathstr,~,~] = fileparts(FileName);
if ~isempty(pathstr) && ~exist(pathstr, 'dir'); mkdir(pathstr);end

Fout = fopen(FileName, 'w'); % use 'a' instead for append

if nargin == 4
    fprintf(Fout, 'IDs');
    fprintf(Fout, '\t%s', ColLabels{:});
    fprintf(Fout, '\n');
    
    FormatS = ['%s' repmat('\t%G', 1, size(Data,2)) '\n'];
    OutBuff = [RowLabels num2cell(Data)]'; clear Data
    fprintf(Fout, FormatS, OutBuff{:,:});
elseif nargin == 2
    %% treat is as Structure with all fields as vercotors of same length.
    Fnames = fieldnames(Data);
    OutBuff = cell(length(Fnames),length(Data.(Fnames{1})));
    FormatS = [];
    GuideFormatS = ['#'];
    fprintf(Fout, ['%s' repmat('\t%s', 1, length(Fnames)-1) '\n'], Fnames{:});
    
    for f = 1: length(Fnames)
        fCol = Data.(Fnames{f});
        if isnumeric(fCol) || islogical(fCol)
            % treat it as number
            OutBuff(f,:) = num2cell(fCol);
            FSchr = 'G';
            GFSchr= '1';
            if isequaln(round(fCol), fCol) || islogical(fCol)
                FSchr = 'd'; 
            end
            
        else
            % Assume it's a string
            OutBuff(f,:) = fCol;
            FSchr = 's';
            GFSchr= 'str';
        end
        FormatS = [FormatS '%' FSchr '\t'];
        GuideFormatS = [GuideFormatS GFSchr '\t'];
    end
    FormatS(end) = 'n';
    GuideFormatS(end) = 'n';
    
    fprintf(Fout, GuideFormatS);
    fprintf(Fout, FormatS, OutBuff{:,:});
end
fclose(Fout);
disp(['Wrote: ' FileName]);

end
