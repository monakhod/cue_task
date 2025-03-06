function BehDataLoadPriorConfidence()
% Load up each eeg data file, extract behaviour, save as one file
% set flagged trials to NaN here
% also save a version of behData that has pps combined across 3 sessions
% 

% clear; close all;

%% 

outFolder = 'D:\cue_task\analysis\Data\Saves';

% loop through each file again, and load the raw eeg file
load(fullfile(outFolder, 'ExtractEpochsPriorConfidence.mat'));

% [name in loaded file, name for behData]
fNames = {  'respLR', 'respLR';
            'RS', 'RS';
            'corrLR', 'corrLR';
            'initAcc', 'initAcc';
            'initResp', 'initResp';
            'blockNum', 'block';
            'cues', 'cues';
            'certainty', 'certainty';
            };
nF = size(fNames,1);

% get these from dataMat
dNames = {
            'certainty2', 'certainty';
            'reward','reward';
            'pp','pp';
            'block','block';
            'day','day';
            'miniBlock','miniBlock';
            'deltaC','deltaC'
            'isBadTrial','isBadTr';
            'gratingFreq', 'gratingFreq';
            'targetFreq', 'targetFreq';
            'chosenFreq', 'chosenFreq';
            };


% prealloc
dataMat = NaN(fileInfo.nPP, fileInfo.maxTr, 18);
for iF = 1:nF
    behData.(fNames{iF,2}) = NaN(fileInfo.nPP, fileInfo.maxTr);
end


for iPP = 1:fileInfo.nPP
    disp(fileInfo.ppID{iPP});
    
    eegFile = fullfile(fileInfo.behFolder, [fileInfo.ppID{iPP} '_beh.mat']);
    
    if exist(eegFile, 'file')
        fprintf('\n loading %s\n', fileInfo.ppID{iPP});
        d = load(eegFile, 'dataMat','dataMatNames', fNames{:,1});
        
        nTr = length(d.RS);
        for iF = 1:nF
            behData.(fNames{iF,2})(iPP,1:nTr) = d.(fNames{iF,1});
        end
        
        behData.pp(iPP,1:nTr) = repmat(iPP,1,nTr); % ppID
%         behData.pp1(iPP,1:nTr) = str2double(fileInfo.ppID1{iPP}(2:3));
        behData.pp1(iPP,1:nTr) = str2double(fileInfo.ppID{iPP}(regexp(fileInfo.ppID{iPP}, '_')+(-2:-1)));
        dataMat(iPP,1:nTr,:) = permute(d.dataMat,[3,1,2]); % store all other data in matrix
        
        behData.trial(iPP,1:nTr) = 1:nTr; % trial number
%         behData.session(iPP,1:nTr) = repmat(str2double(fileInfo.ppID{iPP}(5)),1,nTr);

    else
        fprintf('\n Not Found: %s\n', fileInfo.ppID{iPP});
        
    end
    
end
    
dataMatNames = d.dataMatNames; 

% extract from dataMat too
for iF = 1:size(dNames,1)
    ind = strcmp(dataMatNames, dNames{iF,1});
    if isempty(ind); error('field %s not found', dNames{iF,1}); end
    behData.(dNames{iF,2}) = dataMat(:,:,ind);
end

%% convert

% certainty/uncertainty [maybe probably certain]
behData.certainty = abs(round(behData.initResp -3.5)); %1-6 -> [3 2 1 1 2 3]

% get confidence in correct option
behData.confInCorr = behData.initResp;
behData.confInCorr(behData.corrLR == 1) = 7 - behData.confInCorr(behData.corrLR==1);


%% swap cue numbers so it goes I:N:V, i.e. linear effect

behData.cues(behData.cues == 3) = -1; % invalid
behData.cues(behData.cues == 1) = 0; % neutral
behData.cues(behData.cues == 2) = 1; % valid

%% also get RT in time from samples
behData.RT = behData.RS/eeg.fs * 1000;

%% align all NaNs

nanTr = isnan(behData.RS) | behData.initResp==0;
fn = fieldnames(behData);

for i = 1:length(fn)
    behData.(fn{i})(nanTr) = NaN;
end

dataMat(repmat(nanTr,1,1,length(dataMatNames))) = NaN;

%% set flagged trials to NaN
% 
% load(fullfile(outFolder, 'FlagArtefactsPriorConfidence.mat'), 'isFlagged'); %[nTr nPP] so rotate it
% 
% % for each field
% for i = 1:length(fn)
%     behData.(fn{i})(isFlagged') = NaN; % rotate it to match size
% end
% 
% dataMat(repmat(isFlagged', 1,1,length(dataMatNames))) = NaN; % this too



%% also combine pps across 3 sessions
% [sess*pp tr] into [pp tr*sess]

nPPUniq = length(unique(fileInfo.ppID1n)); % num unique pp - missing numbers will lead to NaN rows

combine2d = @(x) reshape(permute(groupMeans(x, 1, fileInfo.ppID1n,'dim'),[2,1,3]),nPPUniq,[]); %[pp*sess tr] into [pp tr*sess]

behDataPp = structfun(combine2d, behData, 'Uni',0); %[pp tr*sess]

% switch the names over
behDataSess = behData; % each sess + pp is separate row
behData = behDataPp; % each pp gets 1 row, columns go sess:tr
clear behDataPp;

%% labels for things

labels.initAcc = {'error','correct'};
labels.button = {'1','2','3','4','5','6'};
labels.certainty = {'maybe','probably','certain'};
labels.confInCorr = {'certain incorrect', 'probably incorrect', 'maybe incorrect', 'maybe correct', 'probably correct', 'certain correct'};
labels.confResp = {'certain left', 'probably left', 'maybe left', 'maybe right', 'probably right', 'certain right'};
labels.hand = {'left','right'};

labels.cues = {'invalid','neutral','valid'};

%% make table
% use behData so col() goes through [pp tr sess]
% when adding to table, must use the session-combined data e.g. [14 2880]

behVarNames = {'initAcc','RT','certainty','confInCorr'};

nBehVars = length(behVarNames);
behDataNames = fieldnames(behData); % names of all behData DVs

% also get these covars and IVs
behNames = ['pp','pp1','day','block','deltaC', 'trial','cues','isBadTr',...
    'corrLR','miniBlock', 'gratingFreq','targetFreq','chosenFreq', ...
    behVarNames]; % not DVs
behTab = struct2table(structfun(@col, keepFields(behData, behNames), 'UniformOutput',0));
behNames =  behTab.Properties.VariableNames; % update

% make list of names with 'logistic' appended
logIVs = {'initAcc'};
ivNames = behVarNames;
logisticInds = find(ismember(behVarNames, logIVs));
ivNames(logisticInds) = strcat(ivNames(logisticInds), 'Logistic');

% make copy before nanzscoring
for j = logisticInds
    behTab.(ivNames{j}) = behTab.(behNames{strcmp(behNames, behVarNames{j})});
end

% also get log RT 
behTab.RTLog = log(behTab.RT);

behNames = behTab.Properties.VariableNames; % update

isLogistic = ismember(behNames, ivNames(logisticInds)); % which are [0 1] for logistic regressions

glmeArgs = { {}, {'link','logit','distribution','binomial'} }; % for normal/logistic regression, if behVars are DV

% % can z-score now, as will not be excluding trials later?
% behTab(:,~isLogistic) = varfun(@nanzscore, behTab(:,~isLogistic));


%% save

clear d; % don't save this
save(fullfile(fileInfo.behFolder, 'BehDataLoadPriorConfidence.mat'));
save(fullfile(fileInfo.outFolder, 'BehDataLoadPriorConfidence.mat'));
        
end