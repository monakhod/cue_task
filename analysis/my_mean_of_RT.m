clc
clear all

folder = 'C:\Users\monakhvd\OneDrive - Trinity College Dublin\Data\General';
f = what(folder);
% files = f.mat(cellRegexpi(f.mat, '15.23')>0);
files = f.mat(cellfun(@(x) ~isempty(regexpi(x, 'sorted_t')), f.mat));

labels.initAcc = {'error','correct'};
labels.button = {'1','2','3','4','5','6'};
labels.cue = {'invalid','neutral','valid'}; % in order groupMeans returns
labels.certainty = {'maybe','probably','certain'};
labels.cert_cue = {'maybe invalid','probably invalid','certain invalid', ...
    'maybe neutral','probably neutral','certain neutral', 'maybe valid','probably valid','certain valid'};
labels.confInCorr = {'certain incorrect', 'probably incorrect', 'maybe incorrect', 'maybe correct', 'probably correct', 'certain correct'};
labels.initResp = {'certain left', 'probably left', 'maybe left', 'maybe right', 'probably right', 'certain right'};
labels.hand = {'left','right'};
labels.freq = {'20Hz','25Hz'};

labels.cueInitAcc = {'error invalid', 'error neutral', 'error valid',   ...
    'correct invalid', 'correct neutral', 'correct valid'};
combined_labels ={'invalid error maybe' , 'invalid error probably', 'invalid error certain',...
    'invalid correct maybe', 'invalid correct probably', 'invalid correct certain', ...
    'neutral error maybe', 'neutral error probably', 'neutral error certain', 'neutral correct maybe', 'neutral correct probably',...
    'neutral correct certain', 'valid error maybe', 'valid error probably', 'valid error certain',...
    'valid correct maybe', 'valid correct probably', 'valid correct certain'};



resp_1 = struct();
% Loop through each label and add it as a field to the struct
for i = 1:length(labels.initAcc)
    % Replace spaces with underscores to create a valid field name
    field_name = strrep(labels.initAcc{i}, ' ', '_');
    
    % Add the field to the struct with an empty cell as the value
    resp_1.(field_name) = [];
end

resp_2 = struct();
% Loop through each label and add it as a field to the struct
for i = 1:length(labels.cue)
    % Replace spaces with underscores to create a valid field name
    field_name = strrep(labels.cue{i}, ' ', '_');
    
    % Add the field to the struct with an empty cell as the value
    resp_2.(field_name) = [];
end

resp_3 = struct();
% Loop through each label and add it as a field to the struct
for i = 1:length(labels.certainty)
    % Replace spaces with underscores to create a valid field name
    field_name = strrep(labels.certainty{i}, ' ', '_');
    
    % Add the field to the struct with an empty cell as the value
    resp_3.(field_name) = [];
end

resp_4 = struct();
% Loop through each label and add it as a field to the struct
for i = 1:length(labels.cert_cue)
    % Replace spaces with underscores to create a valid field name
    field_name = strrep(labels.cert_cue{i}, ' ', '_');
    
    % Add the field to the struct with an empty cell as the value
    resp_4.(field_name) = [];
end

resp_5 = struct();
% Loop through each label and add it as a field to the struct
for i = 1:length(labels.cueInitAcc)
    % Replace spaces with underscores to create a valid field name
    field_name = strrep(labels.cueInitAcc{i}, ' ', '_');
    
    % Add the field to the struct with an empty cell as the value
    resp_5.(field_name) = [];
end


% Initialize an empty struct
resp_st = struct();

% Loop through each label and add it as a field to the struct
for i = 1:length(combined_labels)
    % Replace spaces with underscores to create a valid field name
    field_name = strrep(combined_labels{i}, ' ', '_');
    
    % Add the field to the struct with an empty cell as the value
    resp_st.(field_name) = [];
end


for i = 1:length (files)
    fprintf('Loading files: %s\n ', files{i});
    load (fullfile(folder,  files{i}));

    %for resp_error/correct
for i = 1: length (labels.initAcc)
 a = nanmean (RT_p(e_or_c ==(i-1))); %finding erp for his condition + mean
field_name = strrep(labels.initAcc{i}, ' ', '_'); %definding name of folder
    
   if ~isfield(resp_1, field_name)
        resp_1.(field_name) = []; 
    end
    
    % Append the computed mean value as a new row to the field
    resp_1.(field_name) = [resp_1.(field_name); a];
end

%for resp and cue
for i = 1: length (labels.cue)
 a = nanmean (RT_p (cue_m ==(i-2))); %finding erp for his condition + mean
field_name = strrep(labels.cue{i}, ' ', '_'); %definding name of folder
    
   if ~isfield(resp_2, field_name)
        resp_2.(field_name) = []; 
    end
    
    % Append the computed mean value as a new row to the field
    resp_2.(field_name) = [resp_2.(field_name); a];
end

%for resp and cert
for i = 1: length (labels.certainty)
 a = nanmean (RT_p(  cert ==i)); %finding erp for his condition + mean
field_name = strrep(labels.certainty{i}, ' ', '_'); %definding name of folder
    
   if ~isfield(resp_3, field_name)
        resp_3.(field_name) = []; 
    end
    
    % Append the computed mean value as a new row to the field
    resp_3.(field_name) = [resp_3.(field_name); a];
end

%for resp and cue+cert

for i = 1: length (labels.cert_cue)
 a = nanmean (RT_p( cue_and_cert ==i)); %finding erp for his condition + mean
field_name = strrep(labels.cert_cue{i}, ' ', '_'); %definding name of folder
    
   if ~isfield(resp_4, field_name)
        resp_4.(field_name) = []; 
    end
    
    % Append the computed mean value as a new row to the field
    resp_4.(field_name) = [resp_4.(field_name); a];
end


for i = 1: length (labels.cueInitAcc)
 a = nanmean (RT_p( cue_and_e ==i)); %finding erp for his condition + mean
field_name = strrep(labels.cueInitAcc{i}, ' ', '_'); %definding name of folder
    
   if ~isfield(resp_5, field_name)
        resp_5.(field_name) = []; 
    end
    
    % Append the computed mean value as a new row to the field
    resp_5.(field_name) = [resp_5.(field_name); a];
end

for i = 1: length (combined_labels)
 a = nanmean (RT_p ( many_cond ==i)); %finding erp for his condition + mean
field_name = strrep(combined_labels{i}, ' ', '_'); %definding name of folder
    
   if ~isfield(resp_st, field_name)
        resp_st.(field_name) = []; 
    end
    
    % Append the computed mean value as a new row to the field
    resp_st.(field_name) = [resp_st.(field_name); a];
end
end


t = []; n = [];
structNames = {'resp_1', 'resp_2', 'resp_3', 'resp_4', 'resp_5', 'resp_st'};

% Iterate over each structure
for i = 1:length(structNames)
    % Get the structure by dynamic field reference
    currentStruct = eval(structNames{i});
    
    % Get field names of the current structure
    fields = fieldnames(currentStruct);
    
    % Iterate over each field
    for j = 1:length(fields)
        fieldName = fields{j};
        
        % Assuming that the field contains data in the expected format
        % Extract the data for the current field
        data = currentStruct.(fieldName);
        
     
            % Compute the nanmean along the second dimension (assuming it's a 3D array)
         %   a = nanmean(data(:,k-300:k-100,:), 2); 
            % Append the result to t
            t = [t, data]; % a(:)' ensures it's a row vector
       n = [n, {[fieldName, '_rt']}];
    end
end
id_1 = 1:length (files);
id = files;
variableNames = [{'id'}, n];
table_of_means_RT = array2table ([id_1', t],  'VariableNames',variableNames)

folder1 = 'C:\Users\monakhvd\Desktop\cue_task\analysis\Data\means';
save (fullfile(folder1, ['mean_RT']),...
        'id', 'table_of_means_RT')