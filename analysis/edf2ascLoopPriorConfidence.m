function edf2ascLoopPriorConfidence(folders)
% convert all edf files in a folder to asc
% edf2asc won't run on OD for some reason, so copy them to local, run, copy asc back and delete edf
% folders is cell array of folders to check
% doesn't overwrite existing, and only converts those without asc already


%%

% % this is folder where data is stored
dataFolder = 'D:\cue_task\analysis\Data\'; %C:\Users\monakhvd\Desktop\cue_task\analysis\Data\'; % 1-16

tmpFolder = 'D:\cue_task\analysis\Data\tmp'; % will be copied here so edf2asc can run
if ~exist(tmpFolder,'dir'); mkdir(tmpFolder); end % make if not exist

overwrite = 0; % don't overwrite existing asc

for j = 1:length(folders)
    % find all edf files in folder
    f = dir(fullfile(dataFolder, folders{j}));
    fNames = {f.name}';
    
    toKeep = cellRegexpi(fNames, 'P\d\d_\d_\d\.edf')>0; % keep only edf files
    f = f(toKeep);
    
    % convert each
    for i = 1:length(f)
    
        edfFileName = f(i).name; % just stem of name
        ascFileName = [f(i).name(1:end-4) '.asc']; % just stem of name
    
        fprintf('\n%s',f(i).name);
        
        if ~exist( fullfile(f(i).folder, ascFileName), 'file' )
            fprintf('... converting to asc');
            
            % copy edf to tmp
            copyfile(fullfile(f(i).folder, edfFileName), tmpFolder);
    
            % convert to asc
            edf2ascMat(fullfile(tmpFolder, edfFileName), [], overwrite);
    
            % move asc back
            movefile(fullfile(tmpFolder, ascFileName), f(i).folder);
    
            % delete edf from here
            delete(fullfile(tmpFolder, edfFileName));
    
        end
    
    end
    

end