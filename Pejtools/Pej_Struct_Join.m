% This function joins structures based on a one or more Keys

% ByUnion, says if you want to Join by union (true), or by intersection of
% the keys (false/default);

% Example:
% Attributes = Pej_Struct_Join( 'SAMPID', [{Pej_Read_Table(JoinAttributesFile)}, {Pej_Read_Table(SizeFactorsFile)}]);


function [InputArray] = Pej_Struct_Join(KeysFields, InputArray, ByUnion)
if nargin<3
    % Join by intersect, discard non-shared rows
    ByUnion = false;
end

if ~iscell(KeysFields)
    KeysFields = {KeysFields};
end
%% Make the Key format
KeyFormat = '';
for i = length(KeysFields):-1:1
    if isnumeric(InputArray{1}.(KeysFields{i})) || islogical(InputArray{1}.(KeysFields{i}))
        if isequalwithequalnans(InputArray{1}.(KeysFields{i}), round(InputArray{1}.(KeysFields{i})))
            %print as integer
            KeyFormat= ['%d_' KeyFormat ];
        else
            % print as float
            KeyFormat= ['%f_' KeyFormat ];
        end
        KeyType(i) = 1; % number
    else
        KeyFormat= [KeyFormat '%s_'];
        KeyType(i) = 2; % string
    end
end
KeyFormat = [KeyFormat(1:end-1) '\n'];

%% build the key for each structure
for K = 1:length(InputArray)
    clear KeyBuff
    for i = length(KeysFields):-1:1
        if KeyType(i) == 1
            KeyBuff(:,i)  = num2cell(InputArray{K}.(KeysFields{i}));
        else
            KeyBuff(:,i)  = (InputArray{K}.(KeysFields{i}));
        end
    end
    KeyBuff   = KeyBuff';
    InputArray{K}.Pej_Key = regexp(sprintf(KeyFormat, KeyBuff{:}), '\n', 'split')';
    if isempty(InputArray{K}.Pej_Key{end})
        % regexp has added one extra cell at the end
        InputArray{K}.Pej_Key(end)=[];
    end
    
    if length(unique(InputArray{K}.Pej_Key))<length(InputArray{K}.Pej_Key)
        warning('Error: Your set of columns do not make a unique key!')
    end
    
    %% Make a list of all keys
    if ByUnion == false
        %% join by intersect
        % Here I discard all the keys that are not observed in all samples
        if K ==1
            AllKeys = InputArray{K}.Pej_Key;
        else
            AllKeys = intersect(AllKeys, InputArray{K}.Pej_Key);
        end
    else
        %% join by union
        % Here I keep all the keys even if they are not observed in all samples
        if K ==1
            AllKeys = InputArray{K}.Pej_Key;
        else
            AllKeys = union(AllKeys, InputArray{K}.Pej_Key);
        end
    end
end
%% Compare everything to the All Keys and return the lists in a comparable index
JL = length(AllKeys);
for K = 1:length(InputArray)
    [~, ia, ib] = intersect(AllKeys, InputArray{K}.Pej_Key);
    FN = fieldnames(InputArray{K});
    F  = length(FN);
    for f = 1:F
        if isnumeric(InputArray{K}.(FN{f}))
            t.(FN{f}) =  nan(length(AllKeys), size(InputArray{K}.(FN{f}),2));
        else
            t.(FN{f}) = cell(length(AllKeys), size(InputArray{K}.(FN{f}),2));
            [t.(FN{f}){:}] = deal('NotMatched'); % 
        end
        t.(FN{f})(ia,:) = InputArray{K}.(FN{f})(ib,:);
    end
    t.Pej_Key = AllKeys;
    InputArray{K} = t; clear t
end


end
