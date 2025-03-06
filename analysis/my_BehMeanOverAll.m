% Based on PriorsConfidenceQuickAnalysis
% Do quick analysis and plots on a single data file
%
% requires: 
% Packages: crameri, barwitherr
% My funcs: 
%   CountUnique.m, ksdensities.m, col.m, StatsTableFromLMECells.m,
% from matlib:
%   groupMeans.m, sq.m, nancat.m, errorBarPlot.m, catStruct.m, 
%   ensureStructsAssignable.m, struct2workspace.m, nanzscore.m
% 


% %%% make sure to change line 22, 24
% close all; clear all; clc
% 
% sub = {'P01' 'P02' 'P03' 'P04'  'P05' 'P06' 'P07' 'P08' ...
%     'P10'  'P13'  'P14' 'P15' 'P16' 'P17' 'P18'};

t_gen = []; conf_by_cue = []; cert_by_cue = []; RT_for_resp =[]; RT_for_cue = []; mean_RT = []; in_cue = []; in_fr = []; in_h = []; in_acc = []; RT_for_e_or_v = []; RT_for_cue_and_corr = [];
 amount_cue_i = []; amount_cue_n = []; amount_cue_v = []; RT_for_18= []; in_cue_cert = [];  amount_cue_i_e =[]; amount_cue_i_c=[];
 amount_cue_n_e=[]; amount_cue_n_c=[]; amount_cue_v_e=[]; amount_cue_v_c=[];cert_by_cue_e = [];

%% initiating all needed structures and variables
%labels
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

labels.cueInitAcc = {'error invalid','correct invalid', 'error neutral','correct neutral', 'error valid',   ...
  'correct valid'};


colours.certainty = [ 0.4974,0.7088,0.9559; 0.5098, 0.3608, 0.7882; 0.7804, 0.1412, 0.4000];

colours.cue_and_corr = [0.925, 0.6625, 0.549; .85, .325, .098;  0.5, 0.724, 0.8705; 0, .447, .741;  0.733, 0.837, 0.594; .466, .674, .188];
colours.confInCorr =  [flipud(crameri('roma',6)); ];
colours.initAcc = [0.6350    0.0780    0.1840; 0.4940    0.1840    0.5560; ];
colours.cmap = crameri('vik');
colours.default = [0 .447 .741; .85 .325 .098; .929 .694 .125;
    .494 .184 .556; .466 .674 .188; .301 .745 .933; .635 .078 .184;];
colours.cue = [.85, .325, .098; 0, .447, .741;.466, .674, .188 ];
colours.for_18 = [0.925, 0.6625, 0.549; .85, .325, .098; 0.425, 0.1625, 0.049;0.925, 0.6625, 0.549; .85, .325, .098; 0.425, 0.1625, 0.049;...
    0.5, 0.724, 0.8705; 0, .447, .741; 0, 0.2235, 0.3705; 0.5, 0.724, 0.8705; 0, .447, .741; 0, 0.2235, 0.3705;...
    0.733, 0.837, 0.594; .466, .674, .188; 0.233, 0.337, 0.094; 0.733, 0.837, 0.594; .466, .674, .188; 0.233, 0.337, 0.094]
colours.cue_and_cert = [0.925, 0.6625, 0.549;.85, .325, .098; 0.425, 0.1625, 0.049;0.5, 0.724, 0.8705; 0, .447, .741; 0, 0.2235, 0.3705; 0.733, 0.837, 0.594; .466, .674, .188; 0.233, 0.337, 0.094]
colors = {  [0, .447, .741], %blue
         [.85, .325, .098], %orange
          [.929, .694, .125],% yellow
          [.494, .184, .556], %purple
          [.466, .674, .188], %green
          [.301, .745, .933], %light blue
          [.635, .078, .184],%red
          [.2, .6, .3], %green
          [.8, .2, .5], %pinkish
          [.1, .3, .7] % another blue
[0.5, 0.724, 0.8705] %11  (lighter blue (1))
[0, 0.2235, 0.3705] %12  (darker blue (1))
[0.925, 0.6625, 0.549] %13  (lighter Orange (2))
[0.425, 0.1625, 0.049] %14  (darker Orange (2))
[0.733, 0.837, 0.594] %15  (lighter green (5))
[0.233, 0.337, 0.094] %16  (darker green (5))
[0.7804, 0.1412, 0.4000] %17 black
[0.5098, 0.3608, 0.7882] %18 grey
[0.4974,0.7088,0.9559]}; %19 lighter grey
nTr = [];


%% going through each participant


for i = 1 : length (sub)
%% loading flagged
folder1 = 'G:\cue_task\analysis\Data\General'


f = what(folder1);
%files = f.mat(cellRegexpi(f.mat, 'flagged')>0);
% files = f.mat(cellfun(@(x) ~isempty(regexpi(x, [sub(i), '_whole.*' ])), f.mat));

pattern = [sub{i}, '_whole.*'];
files = f.mat(cellRegexpi(f.mat, pattern)>0);

disp(files);
nF = length(files);
isGood_g = [];

load(fullfile(folder1, files{1}), 'isGood_gen');
reshaped_data = reshape(isGood_gen, 144, 18);
reshaped_data1 = reshaped_data';
    
a = length (find(isGood_gen==1));
   nTr = [nTr; a];
    isGood_g = reshaped_data1;

%% load up the file into the workspace

% option to load mulitple files

folder1 = 'G:\cue_task\analysis\Data\behF';
%folder1 = 'C:\Users\monakhvd\OneDrive - Trinity College Dublin\Data\P10\beh';


f = what(folder1);
% files = f.mat(cellRegexpi(f.mat, '15.23')>0);
pattern = ['.*', sub{i}, '.*2024.*'];
files = f.mat(cellRegexpi(f.mat, pattern)>0);

%files = f.mat(cellRegexpi(f.mat, '2024')>0);
% files = f.mat(cellRegexpi(f.mat, 'CuePrac')>0);

disp(files);
nF = length(files);

clear resp par PTBevent PTBeventT
v = {'resp','PTBevent','PTBeventT','par'};
for iF = 1:nF
    data{iF} = load(fullfile(folder1, files{iF}),v{:});

    if ~isfield(data{iF}.par, 'targetFreq') % get target frequency
        data{iF}.par.targetFreq = 2 - (data{iF}.par.trialLR == data{iF}.par.gratingFreq); %1=20Hz target, 2=25Hz target

        % also re-calc chosen freq as it was wrong
        data{iF}.resp.freq = 2-(data{iF}.resp.LR(1:length(data{iF}.par.gratingFreq)) == data{iF}.par.gratingFreq);
    end


    % if there are diff numbers of each resp field, top them up so they
    % align later
    fn = fieldnames(data{iF}.resp);
    nDV = structfun(@length, data{iF}.resp);
    inds = find(nDV ~= max(nDV));
    for j = 1:length(inds)
        if any(strcmp(fn{inds(j)},{'otherResps'}))
            data{iF}.resp.(fn{inds(j)})(nDV(inds(j))+1:max(nDV)) = {[]};
        else
            data{iF}.resp.(fn{inds(j)})(:,nDV(inds(j))+1:max(nDV)) = NaN;
        end
    end
    data{iF}.resp.sess = repmat(iF, size(data{iF}.resp.time));

    %     % get folder it came from, use as PP
    %     slashInds = regexp(files{iF},'\');
    %     folderBit = files{iF}(1:slashInds(end)-1);
    %     folderInd = find(strcmp(folder, folderBit));
    %     data{iF}.resp.ppSess = repmat(folderInd, size(data{iF}.resp.time));
    %
    %     % keep P0X bit too
    %     data{iF}.resp.pp = repmat(str2double(folderBit(regexp(folderBit,'P\d\d') + [1 2])), size(data{iF}.resp.time)); %1 for P01 etc

%     data{iF}.resp.pp = repmat(str2double(files{iF}(2:3)), size(data{iF}.resp.time)); % P00
    data{iF}.resp.pp = repmat(iF, size(data{iF}.resp.time));

    fn = fieldnames(data{iF}.resp);

    % combine structs, within each field
    if iF==1
        d = data{iF};
    else
        d.resp = catStruct(2, d.resp, data{iF}.resp); % cat struct
        [dest,src] = ensureStructsAssignable(d.par,data{iF}.par); % make them match
        d.par = cat(2, dest, src); % cat struct, not fields
        for iV = 2:length(v)-1 % vectors
            d.(v{iV}) = cat(2, d.(v{iV}), data{iF}.(v{iV})); %
        end

    end


end


struct2workspace(d);

%nTr = length(isGood (i*i:i+2,:)==1)

%% extract into struct

% get these from resp
behNames = {'initAcc', 'perf';
    'RT', 'time';
    'initResp', 'initResp';
    'confidence', 'confidence';
    'certainty', 'certainty';
    'respLR', 'LR';
    'chosenFreq', 'freq';
    'score','reward';
    'sess','sess';
    'block','block';
    'pp','pp';
    'badTrial','badTrial';
    };

behData = struct();
for g = 1:size(behNames,1)
    behData.(behNames{g,1}) = resp.(behNames{g,2})(1,:);
end

% get these from par

parNames = {'leftFreq', 'gratingFreq';
    'corrLR', 'trialLR';
    'cue', 'cues';
    'corrFreq', 'targetFreq';
    };


for g = 1:size(parNames,1)
    behData.(parNames{g,1}) = nancat(2, par.(parNames{g,2}));

    % need to nanpad these ones to match size of resp fields
    if size(behData.(parNames{g,1}),2) < size(behData.RT,2)
        behData.(parNames{g,1}) = [behData.(parNames{g,1}), NaN(size(behData.RT) - [0 size(behData.(parNames{g,1}),2)])];
    end
end


%% swap cue numbers so it goes I:N:V, i.e. linear effect

behData.cue(behData.cue == 3) = -1; % invalid
behData.cue(behData.cue == 1) = 0; % neutral
behData.cue(behData.cue == 2) = 1; % valid


%% calculate others

% init confinr1 (4-6)
behData.initConfInR1 = behData.initResp;
behData.initConfInR1(behData.respLR==1) = 7 - behData.initConfInR1(behData.respLR == 1);

behData.confInCorr = behData.initResp;
behData.confInCorr(behData.corrLR == 1) = 7 - behData.confInCorr(behData.corrLR==1);

% % remove bad trails?
% validTrials = (behData.badTrial==0 & behData.RT>=0.1);
% 
% behNames = fieldnames(behData);
% for i = 1:size(behNames,1)
%     behData.(behNames{i})(:,~validTrials) = [];
% end

% split by pp
% pp is uniqPP, across sess. ppSes is diff for each sess
behData = structfun(@(x) groupMeans(x,2,behData.pp,'dim')', behData, 'UniformOutput',0); %[pp tr]

nPP = size(behData.initAcc,1);

if nPP == 1
    behData = structfun(@(x) repmat(x,2,1), behData,'Uni',0);
    behData.pp(2,:) = 2;
    nPP = 2;
end



%% get stuff

% split by initial accuracy
behDataByInitAcc = structfun(@(x) groupMeans(x, 2, behData.initAcc, 'dim'), behData, 'Uni',0); %[pp acc tr]

behDataByCert = structfun(@(x) groupMeans(x, 2, behData.certainty, 'dim'), behData, 'Uni',0); %[pp initConf(1-3) tr]
behDataByCue = structfun(@(x) groupMeans(x, 2, behData.cue, 'dim'), behData, 'Uni',0); %[pp confInR1 tr]

behDataByConfInCorr = structfun(@(x) groupMeans(x, 2, behData.confInCorr, 'dim'), behData, 'Uni',0); %[pp initConfInCorr(1-6) tr]




%% make Nan everywhere where it's isGood == 0

a = fieldnames(behData);
for j = 1:length (a)
    b = a{j};
    behData.(a{j})(isGood_g ==0) = NaN;
end

%% basic stats

%basic stats
fprintf('\nmean init acc: %.2f %.2f', nanmean(nanmean(behData.initAcc,2),1)); in_acc = [in_acc ;nanmean(nanmean(behData.initAcc,2),1)];
fprintf('\nmean init acc by hand: %.2f %.2f', nanmean(groupMeans(behData.initAcc,2,behData.respLR)));
in_h = [in_h; nanmean(groupMeans(behData.initAcc,2,behData.respLR))];
fprintf('\nmean init acc by grating freq (20 Hz, 25 Hz): %.2f %.2f', nanmean(groupMeans(behData.initAcc, 2, behData.corrFreq)));
in_fr = [in_fr;  nanmean(groupMeans(behData.initAcc, 2, behData.corrFreq))];
fprintf('\nmean init acc by cue: %.2f %.2f %.2f', nanmean(nanmean(behDataByCue.initAcc,3),1));
in_cue = [in_cue; nanmean(nanmean(behDataByCue.initAcc,3),1)];


in_cue_cert = [in_cue_cert; nanmean(behData.initAcc(behData.certainty == 1 & behData.cue == -1)), nanmean(behData.initAcc(behData.certainty == 2 & behData.cue == -1)), nanmean(behData.initAcc(behData.certainty == 3 & behData.cue == -1)),...
    nanmean(behData.initAcc(behData.certainty == 1 & behData.cue == 0)), nanmean(behData.initAcc(behData.certainty == 2 & behData.cue == 0)), nanmean(behData.initAcc(behData.certainty == 3 & behData.cue == 0)),...
    nanmean(behData.initAcc(behData.certainty == 1 & behData.cue == 1)), nanmean(behData.initAcc(behData.certainty == 2 & behData.cue == 1)), nanmean(behData.initAcc(behData.certainty == 3 & behData.cue == 1))]



fprintf('\nmean amount of each type of answer(maybe, prop, cert) (per cent): %.2f %.2f %.2f', length(find(behData.certainty == 1))/nTr(i), length(find(behData.certainty == 2))/nTr(i), length(find(behData.certainty == 3))/nTr(i));
fprintf('\nmean amount of each type of answer(maybe, prop, cert): %.2f %.2f %.2f', length(find(behData.certainty == 1)), length(find(behData.certainty == 2)), length(find(behData.certainty == 1)));


fprintf('\nmean amount of each type of answer for invalid: %.2f %.2f %.2f', length(find(behData.certainty == 1 & behData.cue == -1))/length(find( behData.cue == -1)), length(find(behData.certainty == 2 & behData.cue == -1))/length(find( behData.cue == -1)), length(find(behData.certainty == 3 & behData.cue == -1))/length(find( behData.cue == -1)));
amount_cue_i = [amount_cue_i; length(find(behData.certainty == 1 & behData.cue == -1))/length(find( behData.cue == -1)), length(find(behData.certainty == 2 & behData.cue == -1))/length(find( behData.cue == -1)), length(find(behData.certainty == 3 & behData.cue == -1))/length(find( behData.cue == -1))];
fprintf('\nmean amount of each type of answer for neutral: %.2f %.2f %.2f', length(find(behData.certainty == 1 & behData.cue == 0))/length(find( behData.cue == 0)), length(find(behData.certainty == 2 & behData.cue == 0))/length(find( behData.cue == 0)), length(find(behData.certainty == 3 & behData.cue == 0))/length(find( behData.cue == 0)));
amount_cue_n = [amount_cue_n; length(find(behData.certainty == 1 & behData.cue == 0))/length(find( behData.cue == 0)), length(find(behData.certainty == 2 & behData.cue == 0))/length(find( behData.cue == 0)), length(find(behData.certainty == 3 & behData.cue == 0))/length(find( behData.cue == 0))];
fprintf('\nmean amount of each type of answer for invalid: %.2f %.2f %.2f\n', length(find(behData.certainty == 1 & behData.cue == 1))/length(find( behData.cue == 1)), length(find(behData.certainty == 2 & behData.cue == 1))/length(find( behData.cue == 1)), length(find(behData.certainty == 3 & behData.cue == 1))/length(find( behData.cue == 1)));
amount_cue_v = [amount_cue_v; length(find(behData.certainty == 1 & behData.cue == 1))/length(find( behData.cue == 1)), length(find(behData.certainty == 2 & behData.cue == 1))/length(find( behData.cue == 1)), length(find(behData.certainty == 3 & behData.cue == 1))/length(find( behData.cue == 1))];

amount_cue_i_e = [amount_cue_i_e; length(find(behData.certainty == 1 & behData.cue == -1 & behData.initAcc == 0))/length(find( behData.cue == -1 & behData.initAcc == 0)), ...
    length(find(behData.certainty == 2 & behData.cue == -1 & behData.initAcc == 0))/length(find( behData.cue == -1 & behData.initAcc == 0)), ...
    length(find(behData.certainty == 3 & behData.cue == -1 & behData.initAcc == 0))/length(find( behData.cue == -1 & behData.initAcc == 0))];
amount_cue_i_c = [amount_cue_i_c; length(find(behData.certainty == 1 & behData.cue == -1 & behData.initAcc == 1))/length(find( behData.cue == -1 & behData.initAcc == 1)), ...
    length(find(behData.certainty == 2 & behData.cue == -1 & behData.initAcc == 1))/length(find( behData.cue == -1 & behData.initAcc == 1)), ...
    length(find(behData.certainty == 3 & behData.cue == -1 & behData.initAcc == 1))/length(find( behData.cue == -1 & behData.initAcc == 1))];
amount_cue_n_e = [amount_cue_n_e; length(find(behData.certainty == 1 & behData.cue == 0 & behData.initAcc == 0))/length(find( behData.cue == 0 & behData.initAcc == 0)), ...
    length(find(behData.certainty == 2 & behData.cue == 0 & behData.initAcc == 0))/length(find( behData.cue == 0 & behData.initAcc == 0)), ...
    length(find(behData.certainty == 3 & behData.cue == 0 & behData.initAcc == 0))/length(find( behData.cue == 0 & behData.initAcc == 0))];
amount_cue_n_c = [amount_cue_n_c; length(find(behData.certainty == 1 & behData.cue == 0 & behData.initAcc == 1))/length(find( behData.cue == 0 & behData.initAcc == 1)), ...
    length(find(behData.certainty == 2 & behData.cue == 0 & behData.initAcc == 1))/length(find( behData.cue == 0 & behData.initAcc == 1)), ...
    length(find(behData.certainty == 3 & behData.cue == 0 & behData.initAcc == 1))/length(find( behData.cue == 0 & behData.initAcc == 1))];
amount_cue_v_e = [amount_cue_v_e; length(find(behData.certainty == 1 & behData.cue == 1 & behData.initAcc == 0))/length(find( behData.cue == 1 & behData.initAcc == 0)), ...
    length(find(behData.certainty == 2 & behData.cue == 1 & behData.initAcc == 0))/length(find( behData.cue == 1 & behData.initAcc == 0)), ...
    length(find(behData.certainty == 3 & behData.cue == 1 & behData.initAcc == 0))/length(find( behData.cue == 1 & behData.initAcc == 0))];
amount_cue_v_c = [amount_cue_v_c; length(find(behData.certainty == 1 & behData.cue == 1 & behData.initAcc == 1))/length(find( behData.cue == 1 & behData.initAcc == 1)), ...
    length(find(behData.certainty == 2 & behData.cue == 1 & behData.initAcc == 1))/length(find( behData.cue == 1 & behData.initAcc == 1)), ...
    length(find(behData.certainty == 3 & behData.cue == 1 & behData.initAcc == 1))/length(find( behData.cue == 1 & behData.initAcc == 1))];

fprintf('\nmean RT: %.2f', nanmean(nanmean (behData.RT)));
mean_RT = [mean_RT; nanmean(nanmean (behData.RT))];
fprintf('\nmean RT for each cue (inv, n,v): %.2f %.2f %.2f', nanmean(behData.RT(behData.cue == -1)), nanmean(behData.RT(behData.cue == 0)), nanmean(behData.RT(behData.cue == 1)));
RT_for_cue = [RT_for_cue; nanmean(behData.RT(behData.cue == -1)), nanmean(behData.RT(behData.cue == 0)), nanmean(behData.RT(behData.cue == 1))];
fprintf('\nnmean RT for each response type (maybe, prop,cert): %.2f %.2f %.2f\n', nanmean(behData.RT(behData.certainty == 1)), nanmean(behData.RT(behData.certainty == 2)), nanmean(behData.RT(behData.certainty == 3)));
RT_for_resp = [RT_for_resp; nanmean(behData.RT(behData.certainty == 1)), nanmean(behData.RT(behData.certainty == 2)), nanmean(behData.RT(behData.certainty == 3))];
RT_for_e_or_v =  [RT_for_e_or_v; nanmean(behData.RT(behData.initAcc == 0)), nanmean(behData.RT(behData.initAcc == 1))];
RT_for_cue_and_corr = [RT_for_cue_and_corr; nanmean(behData.RT(behData.cue == -1 & behData.initAcc ==0)), nanmean(behData.RT(behData.cue == -1 & behData.initAcc ==1)),...
    nanmean(behData.RT(behData.cue == 0 & behData.initAcc ==0)), nanmean(behData.RT(behData.cue == 0 & behData.initAcc ==1)),...
    nanmean(behData.RT(behData.cue == 1 & behData.initAcc ==0)), nanmean(behData.RT(behData.cue == 1 & behData.initAcc ==1))];
RT_for_18= [RT_for_18; nanmean(behData.RT(behData.cue == -1 & behData.initAcc ==0 & behData.certainty ==1)), nanmean(behData.RT(behData.cue == -1 & behData.initAcc ==0 & behData.certainty ==2)), nanmean(behData.RT(behData.cue == -1 & behData.initAcc ==0 & behData.certainty ==3)),...
    nanmean(behData.RT(behData.cue == -1 & behData.initAcc ==1  & behData.certainty ==1)), nanmean(behData.RT(behData.cue == -1 & behData.initAcc ==1  & behData.certainty ==2)), nanmean(behData.RT(behData.cue == -1 & behData.initAcc ==1  & behData.certainty ==3)),...
    nanmean(behData.RT(behData.cue == 0 & behData.initAcc ==0 & behData.certainty ==1)), nanmean(behData.RT(behData.cue == 0 & behData.initAcc ==0 & behData.certainty ==2)), nanmean(behData.RT(behData.cue == 0 & behData.initAcc ==0 & behData.certainty ==3)),...
    nanmean(behData.RT(behData.cue == 0 & behData.initAcc ==1 & behData.certainty ==1)), nanmean(behData.RT(behData.cue == 0 & behData.initAcc ==1 & behData.certainty ==2)), nanmean(behData.RT(behData.cue == 0 & behData.initAcc ==1 & behData.certainty ==3)),...
    nanmean(behData.RT(behData.cue == 1 & behData.initAcc ==0 & behData.certainty ==1)), nanmean(behData.RT(behData.cue == 1 & behData.initAcc ==0 & behData.certainty ==2)), nanmean(behData.RT(behData.cue == 1 & behData.initAcc ==0 & behData.certainty ==3)),...
    nanmean(behData.RT(behData.cue == 1 & behData.initAcc ==1 & behData.certainty ==1)), nanmean(behData.RT(behData.cue == 1 & behData.initAcc ==1 & behData.certainty ==2)),nanmean(behData.RT(behData.cue == 1 & behData.initAcc ==1 & behData.certainty ==3))];


fprintf('\nmean certainty for each cue type (inv, n,v): %.2f %.2f %.2f',nanmean(behData.certainty(behData.cue == -1)), nanmean(behData.certainty(behData.cue == 0)), nanmean(behData.certainty(behData.cue == 1)));
cert_by_cue = [cert_by_cue; nanmean(behData.certainty(behData.cue == -1)), nanmean(behData.certainty(behData.cue == 0)), nanmean(behData.certainty(behData.cue == 1))];
fprintf('\nmean certainty for incorrect responses for each cue type (inv, n,v): %.2f %.2f %.2f',nanmean(behData.certainty(behData.cue == -1 & behData.initAcc == 0)), nanmean(behData.certainty(behData.cue == 0 & behData.initAcc == 0)), nanmean(behData.certainty(behData.cue == 1 & behData.initAcc == 0)));
fprintf('\nmean certainty for correct responses for each cue type (inv, n,v): %.2f %.2f %.2f\n',nanmean(behData.certainty(behData.cue == -1 & behData.initAcc == 1)), nanmean(behData.certainty(behData.cue == 0 & behData.initAcc == 1)), nanmean(behData.certainty(behData.cue == 1 & behData.initAcc == 1)));

cert_by_cue_e = [cert_by_cue_e; nanmean(behData.certainty(behData.cue == -1 & behData.initAcc == 0)), nanmean(behData.certainty(behData.cue == 0 & behData.initAcc == 0)), nanmean(behData.certainty(behData.cue == 1 & behData.initAcc == 0)), nanmean(behData.certainty(behData.cue == -1 & behData.initAcc == 1)), nanmean(behData.certainty(behData.cue == 0 & behData.initAcc == 1)), nanmean(behData.certainty(behData.cue == 1 & behData.initAcc == 1))]

fprintf('\nmean confidence for each cue type (inv, n,v): %.2f %.2f %.2f',nanmean(behData.confidence(behData.cue == -1)), nanmean(behData.confidence(behData.cue == 0)), nanmean(behData.confidence(behData.cue == 1)));
conf_by_cue = [conf_by_cue; nanmean(behData.confidence(behData.cue == -1)), nanmean(behData.confidence(behData.cue == 0)), nanmean(behData.confidence(behData.cue == 1))];



fprintf('\nmean confidence for incorrect responses for each cue type (inv, n,v): %.2f %.2f %.2f',nanmean(behData.confidence(behData.cue == -1 & behData.initAcc == 0)), nanmean(behData.confidence(behData.cue == 0 & behData.initAcc == 0)), nanmean(behData.confidence(behData.cue == 1 & behData.initAcc == 0)));
fprintf('\nmean confidence for correct responses for each cue type (inv, n,v): %.2f %.2f %.2f\n',nanmean(behData.confidence(behData.cue == -1 & behData.initAcc == 1)), nanmean(behData.confidence(behData.cue == 0 & behData.initAcc == 1)), nanmean(behData.confidence(behData.cue == 1 & behData.initAcc == 1)));



%% making table with all 
n =11;
a = movmean(behData.initAcc,n,2);

%finding unique cues
unique_el = unique (behData.cue); unique_el = unique_el(1:3, :)';
%finding unique  correct/error
unique_acc = unique (behData.initAcc); unique_acc = unique_acc(1:2, :)';
%finding unique certainy
unique_cert = unique (behData.certainty); unique_cert = unique_cert(1:3, :)';


t_table = []; g =1; t_in_cell = []; amount_c = [];

%loop for sorting evething based on cue (inv/n/valid), error/correct and
%certainty (maybe/prob/certain

for i = 1:length (unique_el)
    for j = 1:length (unique_acc)
        for n = 1:length (unique_cert)
            c = sum (sum ( behData.certainty == unique_cert(n) & behData.cue == unique_el(i) & behData.initAcc == unique_acc(j)));
            t = behData.RT( behData.certainty == unique_cert(n) & behData.cue == unique_el(i) & behData.initAcc == unique_acc(j));
            t_m = mean (t);
            amount_c = [amount_c, c];
            t_table = [t_table , t_m];
            g = g+1;
        end
    end
end

t_gen = [t_gen; t_table];





end


combined_labels = cell(length(unique_el) * length(unique_acc)*length(unique_cert), 1);
c = 1;
%making labels for graphs with RT
for l = 1:length(unique_el)
    for j = 1:length(unique_acc)
        for n = 1:length (unique_cert)
        combined_labels{c} = sprintf('%s %s %s', labels.cue{l}, labels.initAcc{j}, labels.certainty{n});
        c = c + 1;
        end
            end
end

%% Diff graphs based on privious table and cell array
x = 1:18; % for x-axes

%plot of sorted mean RT 
figure;
errorBarPlot(t_gen, 'area',1);
hold on
plot(x, polyval(polyfit(x, nanmean (t_gen, 1), 1), x), '--', 'LineWidth', 1, 'color', colors{1});
xline (1:18, '--');
hold off

xticks(1:length(x));
xticklabels(combined_labels);
xtickangle(60); % Rotating x-axis labels for better readability

% Labels and title
xlabel('Category');
ylabel('RT');
title('RT sorted by cue, correct-error and certainty right against left hand');
legend ('mean RT');


%% plotting stuff based on cue and certainty
%labels for graphs
c = 1;
for i = 1:length(unique_el)
        for n = 1:length (unique_cert)
        combined_labels_2{c} = sprintf('%s %s %s', labels.cue{i}, labels.certainty{n});
        c = c + 1;
        end
end 

%% doing Friedman test
g = t_gen (:, [1:3 7:9 13:15]);
k = t_gen (:, [4:6 10:12 16:18]);

s_p = [];
for i = 1:9
    condition1 = g(:, i);
    condition2 = k(:, i);
    
    % Perform Mann-Whitney U test
    p = ranksum(condition1, condition2);
    
    % Determine if the result is significant
    significant_points = p < 0.05;
    s_p = [s_p, significant_points];
end

% Find indices of significant points
indices = find(s_p == 1);

% Display the indices
disp('Indices of significant points:');
disp(indices);

% errorbar 
y = 1:9;
figure
errorBarPlot(g,  'area', 1, 'alpha', 0.2, 'color', colors{2});
hold on
errorBarPlot(k, 'area',1, 'alpha', 0.2, 'color', colors{1});

 plot(y, polyval(polyfit(y, nanmean (g, 1), 1), y), '--', 'LineWidth',1, 'color', colors{2});
 plot(y, polyval(polyfit(y, nanmean (k, 1), 1), y), '--', 'LineWidth', 1, 'color', colors{1});


xline (1:9, '--');
hold off

xticks(1: (length(y)));
xticklabels(combined_labels_2);
xtickangle(60); % Rotating x-axis labels for better readability

% Labels and title
xlabel('Category');
ylabel('RT');
title('RT sorted by cue and certainty');
legend ('', 'RT for error',  '', 'RT for correct', 'trend for error', 'trend for correct');


%% plot for mean RT for error and correct responce
%error +correct

figure
errorBarPlot(g (:,1:3),  'area', 1, 'alpha', 0.2, 'color', colors{2}, 'plotArgs', {'LineStyle', '--'});
hold on
errorBarPlot(k (:, 1:3),  'area', 1, 'alpha', 0.2, 'color', colors{2});

errorBarPlot(g (:,4:6),  'area', 1, 'alpha', 0.2, 'color', colors{1},'plotArgs', {'LineStyle', '--'});
errorBarPlot(k (:,4:6),  'area', 1, 'alpha', 0.2, 'color', colors{1});
errorBarPlot(g (:,7:9),  'area', 1, 'alpha', 0.2, 'color', colors{5}, 'plotArgs', {'LineStyle', '--'});
errorBarPlot(k (:,7:9),  'area', 1, 'alpha', 0.2, 'color', colors{5});




xticks(1: 3);
xticklabels(labels.certainty);
xtickangle(60);

xlabel('Category');
ylabel('RT');
title('RT sorted by cue ');
legend ('', ' RT for error + invalid','', ' RT for correct + invalid', '', 'RT for error + neutral','', 'RT for correct + neutral', '', 'RT for error + valid',  '', 'RT for correct + valid');


%error
figure
errorBarPlot(g (:,1:3),  'area', 1, 'alpha', 0.2, 'color', colors{2});
hold on
errorBarPlot(g (:,4:6),  'area', 1, 'alpha', 0.2, 'color', colors{1});
errorBarPlot(g (:,7:9),  'area', 1, 'alpha', 0.2, 'color', colors{5});

xticks(1: 3);
xticklabels(labels.certainty);
xtickangle(60);

xlabel('Category');
ylabel('RT');
title('RT sorted by cue for error trails only');
legend ('', ' RT for error + invalid', '', 'RT for error + neutral', '', 'RT for error + valid');


%correct
figure
errorBarPlot(k (:, 1:3),  'area', 1, 'alpha', 0.2, 'color', colors{2});
hold on
errorBarPlot(k (:,4:6),  'area', 1, 'alpha', 0.2, 'color', colors{1});
errorBarPlot(k (:,7:9),  'area', 1, 'alpha', 0.2, 'color', colors{5});

xticks(1: 3);
xticklabels(labels.certainty);
xtickangle(60);

xlabel('Category');
ylabel('RT');
title('RT sorted by cue for correct trails only');
legend ('', ' RT for correct + invalid', '', 'RT for correct + neutral', '', 'RT for correct + valid');


%error
figure
errorBarPlot(g (:,1:3:9),  'area', 1, 'alpha', 0.2, 'color', colors{19});
hold on
errorBarPlot(g (:,2:3:9),  'area', 1, 'alpha', 0.2, 'color', colors{18});
errorBarPlot(g (:,3:3:9),  'area', 1, 'alpha', 0.2, 'color', colors{17});

xticks(1: 3);
xticklabels(labels.cue);
xtickangle(60);

xlabel('Category');
ylabel('RT');
title('RT sorted by cue for error trails only');
legend ('', ' RT for error + maybe', '', 'RT for error + probably', '', 'RT for error + certain');


%correct
figure
errorBarPlot(k (:,1:3:9),  'area', 1, 'alpha', 0.2, 'color', colors{19});
hold on
errorBarPlot(k (:,2:3:9),  'area', 1, 'alpha', 0.2, 'color', colors{18});
errorBarPlot(k (:,3:3:9),  'area', 1, 'alpha', 0.2, 'color', colors{17});

xticks(1: 3);
xticklabels(labels.cue);
xtickangle(60);

xlabel('Category');
ylabel('RT');
title('RT sorted by cue for correct trails only');
legend ('', ' RT for correct + maybe', '', 'RT for correct + probably', '', 'RT for correct + certain');


%% simple graphs

fprintf('\nmean init acc: %.2f %.2f', nanmean(in_acc)); 
fprintf('\nmean init acc by hand: %.2f %.2f', nanmean(in_h));
fprintf('\nmean init acc by grating freq (20 Hz, 25 Hz): %.2f %.2f', nanmean(in_fr)); 
fprintf('\nmean init acc by cue: %.2f %.2f %.2f', nanmean(in_cue,1));

fprintf('\namount of answers for invalid cue: %.2f %.2f %.2f', nanmean(amount_cue_i,1));
fprintf('\namount of answers for neutral cue: %.2f %.2f %.2f', nanmean(amount_cue_n,1));
fprintf('\namount of answers for valid cue: %.2f %.2f %.2f', nanmean(amount_cue_v,1));


fprintf('\namount of answers for invalid cue +error: %.2f %.2f %.2f', nanmean(amount_cue_i_e,1));
fprintf('\namount of answers for invalid cue +correct: %.2f %.2f %.2f', nanmean(amount_cue_i_c,1));
fprintf('\namount of answers for neutral cue +error: %.2f %.2f %.2f', nanmean(amount_cue_n_e,1));
fprintf('\namount of answers for neutral cue +correct: %.2f %.2f %.2f', nanmean(amount_cue_n_c,1));
fprintf('\namount of answers for valid cue+error: %.2f %.2f %.2f', nanmean(amount_cue_v_e,1));
fprintf('\namount of answers for valid cue+correct: %.2f %.2f %.2f', nanmean(amount_cue_v_c,1));

fprintf('\nmean RT: %.2f', nanmean(mean_RT));
fprintf('\nmean RT for each cue (inv, n,v): %.2f %.2f %.2f', nanmean(RT_for_cue));
fprintf('\nnmean RT for each response type (maybe, prop,cert): %.2f %.2f %.2f\n', nanmean(RT_for_resp));

fprintf('\nmean certainty for each cue type (inv, n,v): %.2f %.2f %.2f',nanmean(cert_by_cue));
fprintf('\nmean confidence for each cue type (inv, n,v): %.2f %.2f %.2f',nanmean(conf_by_cue));


%% simple graphs
n = 11; % smoothing window
figure();
t = tiledlayout('flow');


% mean acc
ax = nexttile(t);
h = errorBarPlot(in_acc, 'type','bar');
ylim([.5 1]);
ylabel('mean acc');


% resp hand
ax = nexttile(t);
h = errorBarPlot(in_h, 'type','bar');
xticks(1:2);
xticklabels(labels.hand);
ylim([.5 1]);
ylabel('mean acc');
% legend(labels.hand,'Location','North');


% correct freq
ax = nexttile(t);
h = errorBarPlot(in_fr, 'type','bar');
xticks(1:2);
xticklabels(labels.freq);
ylim([.5 1]);
ylabel('mean acc');
 legend('20 Hz',  '25Hz', '','Location','North');

 %cue + accuracy
ax = nexttile(t);
h = errorBarPlot(in_cue, 'type','bar');
h.FaceColor = 'flat';
h.CData = colours.cue;
xticklabels(labels.cue)
xlabel('cue type');
ylabel('accuracy');


 % certainty + cue
ax = nexttile(t);
h = errorBarPlot(cert_by_cue, 'type','bar');
h.FaceColor = 'flat';
h.CData = colours.cue;
xticklabels(labels.cue)
xlabel('cue type');
ylabel('certainty');
ylim([2 2.7]);

% Add value labels on top of the bars
yData = nanmean(cert_by_cue);  % Get the mean data used for the bars
xData = 1:length(yData);       % X locations of the bars
text(xData, yData, num2str(yData', '%.2f'), ...  % Format numbers to 2 decimal places
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 10);

 % certainty + cue + error
ax = nexttile(t);
h = errorBarPlot(cert_by_cue_e, 'type','bar');
h.FaceColor = 'flat';
h.CData = [colours.cue; colours.cue];
xticklabels(labels.cueInitAcc([1, 3, 5, 2, 4, 6]))
xlabel('cue type');
ylabel('certainty');
ylim([2 2.7]);

% Add value labels on top of the bars
yData = nanmean(cert_by_cue_e);  % Get the mean data used for the bars
xData = 1:length(yData);       % X locations of the bars
text(xData, yData, num2str(yData', '%.2f'), ...  % Format numbers to 2 decimal places
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 10);

 % confidence + cue
ax = nexttile(t);
h = errorBarPlot(conf_by_cue, 'type','bar');
h.FaceColor = 'flat';
h.CData = colours.cue;
xticklabels(labels.cue)
xlabel('cue type');
ylabel('confidence');
ylim([5.3 5.6]);

% Add value labels on top of the bars
yData = nanmean(conf_by_cue);  % Get the mean data used for the bars
xData = 1:length(yData);       % X locations of the bars
text(xData, yData, num2str(yData', '%.2f'), ...  % Format numbers to 2 decimal places
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 10);

ax = nexttile(t);
h = errorBarPlot (RT_for_cue, 'type','bar');
h.FaceColor = 'flat';
h.CData = colours.cue;
xticklabels(labels.cue)
xlabel('cue type');
ylabel('mean RT');
ylim([0.7 1]);

% Add value labels on top of the bars
yData = nanmean(RT_for_cue);  % Get the mean data used for the bars
xData = 1:length(yData);       % X locations of the bars
text(xData, yData, num2str(yData', '%.2f'), ...  % Format numbers to 2 decimal places
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 10);


 % confidence + cue
ax = nexttile(t);
h = errorBarPlot(RT_for_resp, 'type','bar');
h.FaceColor = 'flat';
h.CData = colours.certainty;
xticklabels(labels.certainty)
xlabel('response type');
ylabel('mean RT');
ylim([0.7 1.1]);

% Add value labels on top of the bars
yData = nanmean(RT_for_resp);  % Get the mean data used for the bars
xData = 1:length(yData);       % X locations of the bars
text(xData, yData, num2str(yData', '%.2f'), ...  % Format numbers to 2 decimal places
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 10);


% mean RT for error and correct
ax = nexttile(t);
h = errorBarPlot(RT_for_e_or_v, 'type','bar');
h.FaceColor = 'flat';
h.CData = colours.certainty;
xticklabels(labels.initAcc)
%xlabel('');
ylabel('mean RT');
ylim([0 1.3]);

% Add value labels on top of the bars
yData = nanmean(RT_for_e_or_v);  % Get the mean data used for the bars
xData = 1:length(yData);       % X locations of the bars
text(xData, yData, num2str(yData', '%.2f'), ...  % Format numbers to 2 decimal places
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 10);

% mean RT from each cue for error and correct
labels.cue_and_corr = {'error invalid','correct invalid',  'error neutral', 'correct neutral','error valid',   ...
      'correct valid'};

ax = nexttile(t);
h = errorBarPlot(RT_for_cue_and_corr, 'type','bar');
h.FaceColor = 'flat';
h.CData = colours.cue_and_corr ;
xticklabels(labels.cue_and_corr)
%xlabel('');
ylabel('mean RT');
ylim([0.7 1.05]);

% Add value labels on top of the bars
yData = nanmean(RT_for_cue_and_corr);  % Get the mean data used for the bars
xData = 1:length(yData);       % X locations of the bars
text(xData, yData, num2str(yData', '%.2f'), ...  % Format numbers to 2 decimal places
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 10);


%%% new plot of rt by 18 vars

labels.for_18= {'error invalid maybe', 'error invalid probably','error invalid certain',...
    'correct invalid maybe', 'correct invalid probably', 'correct invalid certain',...
    'error neutral maybe', 'error neutral probably', 'error neutral certain',...
    'correct neutral maybe', 'correct neutral probably', 'correct neutral certain',...
    'error valid maybe',  'error valid probably' 'error valid certain'...
      'correct valid maybe', 'correct valid probably', 'correct valid certain'};

ax = nexttile(t);
h = errorBarPlot(RT_for_18, 'type','bar');
h.FaceColor = 'flat';
h.CData = colours.for_18 ;
xticklabels(labels.for_18)
%xlabel('');
ylabel('mean RT');
ylim([0.7 1.1]);

% Add value labels on top of the bars
yData = nanmean(RT_for_18);  % Get the mean data used for the bars
xData = 1:length(yData);       % X locations of the bars
text(xData, yData, num2str(yData', '%.2f'), ...  % Format numbers to 2 decimal places
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 10);

%accuracy for each cue and certanty type
labels.cue_and_cert = {'maybe invalid','probably invalid', 'certain invalid',  'maybe neutral', 'probably neutral', 'certain neutral','maybe valid',   ...
      'probably valid', 'certain valid'};

ax = nexttile(t);
h = errorBarPlot(in_cue_cert, 'type','bar');
h.FaceColor = 'flat';
h.CData = colours.cue_and_cert ;
xticklabels(labels.cue_and_cert)
%xlabel('cue+cert');
ylabel('accuracy');
ylim([0.4 1]);

% Add value labels on top of the bars
yData = nanmean(in_cue_cert);  % Get the mean data used for the bars
xData = 1:length(yData);       % X locations of the bars
text(xData, yData, num2str(yData', '%.2f'), ...  % Format numbers to 2 decimal places
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 10);

%plot for each person
ax = nexttile(t);
h = plot( RT_for_resp', 'o-'); 
   xlabel('response type');
ylabel('mean RT');
xticklabels(labels.certainty)
%legend(sub);

%plot for each person of accuracy for cue
ax = nexttile(t);
h = plot( in_cue', 'o-'); 
xticklabels(labels.cue)
xlabel('cue type');
ylabel('accuracy');
%% additional figures

figure
h1 = errorBarPlot (amount_cue_i_e);
h1.Color = colours.cue(1, :);
h1.LineWidth = 2;
h1.LineStyle = "--";
hold on
h2 =  errorBarPlot (amount_cue_i_c)
h2.Color = colours.cue(1, :)
h2.LineWidth = 2;
h3 =  errorBarPlot (amount_cue_n_e)
h3.Color = colours.cue(2, :);
h3.LineWidth = 2;
h3.LineStyle = "--";
h4 =  errorBarPlot (amount_cue_n_c);
h4.Color = colours.cue(2, :);
h4.LineWidth = 2;
h5 =  errorBarPlot (amount_cue_v_e);
h5.Color = colours.cue(3, :);
h5.LineWidth = 2;
h5.LineStyle = "--";
h6 =  errorBarPlot (amount_cue_v_c);
h6.Color = colours.cue(3, :);
h6.LineWidth = 2;
set(gca, 'XTick', 1:length(labels.cue)); % Assuming x-ticks are in position 1, 2, 3
set(gca, 'XTickLabel', labels.certainty);
ylabel('% of answers ');
xlabel('certainty');
legend (labels.cueInitAcc)
title('Percentage of Responses by Answer Type Across Different Cues');


labels.for_18= {'error invalid maybe', 'error invalid probably','error invalid certain',...
    'correct invalid maybe', 'correct invalid probably', 'correct invalid certain',...
    'error neutral maybe', 'error neutral probably', 'error neutral certain',...
    'correct neutral maybe', 'correct neutral probably', 'correct neutral certain',...
    'error valid maybe',  'error valid probably' 'error valid certain'...
      'correct valid maybe', 'correct valid probably', 'correct valid certain'};

figure;
h = errorBarPlot(RT_for_18, 'type','bar');
h.FaceColor = 'flat';
h.CData = colours.for_18 ;
xticklabels(labels.for_18)
%xlabel('');
ylabel('mean RT');
ylim([0.7 1.1]);
set(gca, 'XTick', 1:length(labels.for_18)); % Assuming x-ticks are in position 1, 2, 3
set(gca, 'XTickLabel', labels.for_18);
% Add value labels on top of the bars
yData = nanmean(RT_for_18);  % Get the mean data used for the bars
xData = 1:length(yData);       % X locations of the bars
text(xData, yData, num2str(yData', '%.2f'), ...  % Format numbers to 2 decimal places
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 10);

%cert+cue+acc

figure;
h = errorBarPlot(cert_by_cue_e, 'type','bar');
h.FaceColor = 'flat';
h.CData = [colours.cue; colours.cue];
xticklabels(labels.cueInitAcc([1, 3, 5, 2, 4, 6]))
xlabel('cue type');
ylabel('certainty');
ylim([2 2.7]);

% Add value labels on top of the bars
yData = nanmean(cert_by_cue_e);  % Get the mean data used for the bars
xData = 1:length(yData);       % X locations of the bars
text(xData, yData, num2str(yData', '%.2f'), ...  % Format numbers to 2 decimal places
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 10);
