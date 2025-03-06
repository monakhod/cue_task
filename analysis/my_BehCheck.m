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


%%% make sure to change line 22, 24
%close all; clear all; clc

%% load up the file into the workspace

% option to load mulitple files

% folder1 = 'C:\Users\monakhvd\OneDrive - Trinity College Dublin\Data\P01\beh';
% 
% wh = 'P18';
f = what(folder1);
% files = f.mat(cellRegexpi(f.mat, '15.23')>0);
files = f.mat(cellRegexpi(f.mat, '2024')>0);
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

nTr = sum(resp.badTrial==0)

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
for i = 1:size(behNames,1)
    behData.(behNames{i,1}) = resp.(behNames{i,2})(1,:);
end

% get these from par

parNames = {'leftFreq', 'gratingFreq';
    'corrLR', 'trialLR';
    'cue', 'cues';
    'corrFreq', 'targetFreq';
    };


for i = 1:size(parNames,1)
    behData.(parNames{i,1}) = nancat(2, par.(parNames{i,2}));

    % need to nanpad these ones to match size of resp fields
    if size(behData.(parNames{i,1}),2) < size(behData.RT,2)
        behData.(parNames{i,1}) = [behData.(parNames{i,1}), NaN(size(behData.RT) - [0 size(behData.(parNames{i,1}),2)])];
    end
end

%colors
colors = cell(1, 18);

%making color
% Darker to lighter red
colors{1} = [0.7, 0, 0];    % Dark Red
colors{2} = [1, 0, 0];      % Medium Red
colors{3} = [1, 0.7, 0.7];  % Light Red

% Darker to lighter blue
colors{4} = [0, 0, 0.7];    % Dark Blue
colors{5} = [0, 0, 1];      % Medium Blue
colors{6} = [0.7, 0.7, 1];  % Light Blue

% Darker to lighter orange
colors{7} = [0.8, 0.4, 0];    % Dark Orange
colors{8} = [1, 0.5, 0];      % Medium Orange
colors{9} = [1, 0.8, 0.5];    % Light Orange

% Darker to lighter green
colors{10} = [0, 0.5, 0];      % Dark Green
colors{11} = [0, 1, 0];        % Medium Green
colors{12} = [0.5, 1, 0.5];    % Light Green

% Darker to lighter yellow
colors{13} = [0.8, 0.8, 0];    % Dark Yellow
colors{14} = [1, 1, 0];        % Medium Yellow
colors{15} = [1, 1, 0.8];      % Light Yellow

% Darker to lighter violet
colors{16} = [0.5, 0, 0.5];    % Dark Violet
colors{17} = [0.7, 0.4, 0.7];  % Medium Violet
colors{18} = [1, 0.7, 1];  
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

% % remove bad trials?
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

labels.initAcc = {'error','correct'};
labels.button = {'1','2','3','4','5','6'};
labels.cue = {'invalid','neutral','valid'}; % in order groupMeans returns
labels.certainty = {'maybe','probably','certain'};
labels.confInCorr = {'certain incorrect', 'probably incorrect', 'maybe incorrect', 'maybe correct', 'probably correct', 'certain correct'};
labels.initResp = {'certain left', 'probably left', 'maybe left', 'maybe right', 'probably right', 'certain right'};
labels.hand = {'left','right'};
labels.freq = {'20Hz','25Hz'};


colours.certainty = [0.4660    0.6740    0.1880; 0.3010    0.7450    0.9330; 0    0.4470    0.7410; ];
colours.confInCorr =  [flipud(crameri('roma',6)); ];
colours.initAcc = [0.6350    0.0780    0.1840; 0.4940    0.1840    0.5560; ];
colours.cmap = crameri('vik');
colours.default = [0 .447 .741; .85 .325 .098; .929 .694 .125;
    .494 .184 .556; .466 .674 .188; .301 .745 .933; .635 .078 .184;];
colours.cue = colours.default(1:3,:);

%basic stats
fprintf('\nmean init acc: %.2f %.2f', nanmean(nanmean(behData.initAcc,2),1));
fprintf('\nmean init acc by hand: %.2f %.2f', nanmean(groupMeans(behData.initAcc,2,behData.respLR)));
fprintf('\nmean init acc by grating freq (20 Hz, 25 Hz): %.2f %.2f', nanmean(groupMeans(behData.initAcc, 2, behData.corrFreq)));
fprintf('\nmean init acc by cue: %.2f %.2f %.2f', nanmean(nanmean(behDataByCue.initAcc,3),1));

fprintf('\nmean amount of each type of answer(maybe, prop, cert) (per cent): %.2f %.2f %.2f', sum(sum (sum(behData.certainty == 1)))/nTr, sum(sum (sum(behData.certainty == 2)))/nTr, sum(sum (sum(behData.certainty == 3)))/nTr);
fprintf('\nmean amount of each type of answer(maybe, prop, cert): %.2f %.2f %.2f', sum(sum (sum(behData.certainty == 1))), sum(sum (sum(behData.certainty == 2))), sum(sum (sum(behData.certainty == 3))));
fprintf('\nmean amount of each type of answer for invalid: %.2f %.2f %.2f', sum(sum (sum(behData.certainty == 1 & behData.cue == -1))), sum(sum (sum(behData.certainty == 2 & behData.cue == -1))), sum(sum (sum(behData.certainty == 3 & behData.cue == -1))));
fprintf('\nmean amount of each type of answer for neutral: %.2f %.2f %.2f', sum(sum (sum(behData.certainty == 1 & behData.cue == 0))), sum(sum (sum(behData.certainty == 2 & behData.cue == 0))), sum(sum (sum(behData.certainty == 3 & behData.cue == 0))));
fprintf('\nmean amount of each type of answer for invalid: %.2f %.2f %.2f\n', sum(sum (sum(behData.certainty == 1 & behData.cue == 1))), sum(sum (sum(behData.certainty == 2 & behData.cue == 1))), sum(sum (sum(behData.certainty == 3 & behData.cue == 1))));

fprintf('\nmean RT: %.2f', nanmean(nanmean (behData.RT)));
fprintf('\nmean RT for each cue (inv, n,v): %.2f %.2f %.2f', nanmean(behData.RT(behData.cue == -1)), nanmean(behData.RT(behData.cue == 0)), nanmean(behData.RT(behData.cue == 1)));
fprintf('\nnmean RT for each response type (maybe, prop,cert): %.2f %.2f %.2f\n', nanmean(behData.RT(behData.certainty == 1)), nanmean(behData.RT(behData.certainty == 2)), nanmean(behData.RT(behData.certainty == 3)));

fprintf('\nmean certainty for each cue type (inv, n,v): %.2f %.2f %.2f',nanmean(behData.certainty(behData.cue == -1)), nanmean(behData.certainty(behData.cue == 0)), nanmean(behData.certainty(behData.cue == 1)));
fprintf('\nmean certainty for incorrect responses for each cue type (inv, n,v): %.2f %.2f %.2f',nanmean(behData.certainty(behData.cue == -1 & behData.initAcc == 0)), nanmean(behData.certainty(behData.cue == 0 & behData.initAcc == 0)), nanmean(behData.certainty(behData.cue == 1 & behData.initAcc == 0)));
fprintf('\nmean certainty for correct responses for each cue type (inv, n,v): %.2f %.2f %.2f\n',nanmean(behData.certainty(behData.cue == -1 & behData.initAcc == 1)), nanmean(behData.certainty(behData.cue == 0 & behData.initAcc == 1)), nanmean(behData.certainty(behData.cue == 1 & behData.initAcc == 1)));


fprintf('\nmean confidence for each cue type (inv, n,v): %.2f %.2f %.2f',nanmean(behData.confidence(behData.cue == -1)), nanmean(behData.confidence(behData.cue == 0)), nanmean(behData.confidence(behData.cue == 1)));
fprintf('\nmean confidence for incorrect responses for each cue type (inv, n,v): %.2f %.2f %.2f',nanmean(behData.confidence(behData.cue == -1 & behData.initAcc == 0)), nanmean(behData.confidence(behData.cue == 0 & behData.initAcc == 0)), nanmean(behData.confidence(behData.cue == 1 & behData.initAcc == 0)));
fprintf('\nmean confidence for correct responses for each cue type (inv, n,v): %.2f %.2f %.2f\n',nanmean(behData.confidence(behData.cue == -1 & behData.initAcc == 1)), nanmean(behData.confidence(behData.cue == 0 & behData.initAcc == 1)), nanmean(behData.confidence(behData.cue == 1 & behData.initAcc == 1)));

%making a table, where all bad outcomes with blocks accuracy are checked
y = nanmean(behDataByCue.initAcc, 3);
sorted_table = [];

for i = 1 :length (y)
    if y (i,2) <0.58
     sorted_table = [sorted_table ;  {i},  {y(i,1)},{y(i,2)} ,{y(i,3)}, {'accuracy for neutral less then 60%'}];
    end
if y (i,2) >0.84
     sorted_table = [sorted_table ;  {i},  {y(i,1)},{y(i,2)} ,{y(i,3)}, {'accuracy is too high for nuetral cues'}]; 
end
if y(i,3) - y (i,2) < 0.05
     sorted_table = [sorted_table ;  {i},  {y(i,1)},{y(i,2)} ,{y(i,3)}, {'small diff between valid and neutral'}];
 end
    if y(i,1) > y(i,2)
        sorted_table = [sorted_table ; {i},  {y(i,1)},{y(i,2)} ,{y(i,3)}, {'invalid bigger then neutrals'}];
    end
 if y (i,3) == 1
     sorted_table = [sorted_table ;  {i},  {y(i,1)},{y(i,2)} ,{y(i,3)}, {'100% accuracy for valid'}];
 end
    if y (i,1) <0.15
     sorted_table = [sorted_table ;  {i},  {y(i,1)},{y(i,2)} ,{y(i,3)}, {'small accuracy for invalid cues'}]; 
end

end

if ~isempty(sorted_table)
what_should_be_del = []; %Calculate the percentage of responses relative to the total number of occurrences
what_should_be_del = array2table ([sorted_table], 'VariableNames',{ 'block', 'inv', 'neutral', 'valid', 'variant'}) %table with names for them
%keyboard
end

%accuracy per each block
figure
bar (nanmean(behDataByCue.initAcc,3));
set(gca,'ColorOrder',colours.certainty,'nextplot','replacechildren')
xlabel ("blocks");
ylabel ("accuracy");
title('Accuracy for each cue during different blocks');
legend (labels.cue, 'Location','Best');


a = nanmean(behDataByCue.initAcc,3);
a2 = a(:,2);
figure
bar (a2);
set(gca,'ColorOrder',colours.certainty,'nextplot','replacechildren')
xlabel ("DeltaC");
ylabel ("accuracy");

xlabels = [par.deltaC];
set(gca, 'XTickLabel', xlabels);
for i = 1:numel(a2)
    text(i, a2(i), num2str(a2(i)), ...
        'HorizontalAlignment',  'center', ...
        'VerticalAlignment', 'bottom');
end

title('Accuracy for each neutral cue during different blocks');
legend (labels.cue (2), 'Location','Best');

acc_for_n_cue = []; %Calculate the percentage of responses relative to the total number of occurrences
acc_for_n_cue = array2table ([behData.pp(:,1), a2, xlabels'], 'VariableNames',{ 'block', 'accracy for neutral', 'deltaC'}) %table with names for them


figure
plot (groupMeans(behData.initAcc, 2, behData.corrFreq))
xlabel('DeltaC');
ylabel('accuracy');
xticks(1: 18);
xticklabels([par.deltaC]);
xtickangle(60);
%% simple graphs

%Accuracy for each cue during different blocks
n = 11; % smoothing window
figure();
t = tiledlayout('flow');

% accs
ax = nexttile(t);
h = errorBarPlot(movmean(behData.initAcc,n,2), 'area',1);
xlabel('trial');
ylabel('accuracy');

% rt during block
ax = nexttile(t);
h = errorBarPlot(movmean(behData.RT,n,2), 'area',1);
xlabel('trial');
ylabel('RT (s)');

% cert during block
ax = nexttile(t);
errorBarPlot(movmean(behData.certainty,n,2), 'area',1);
xlabel('trial');
ylabel('certainty');

% score
ax = nexttile(t);
errorBarPlot(movmean(behData.score,n,2), 'area',1);
xlabel('trial');
ylabel('score');

% mean acc
ax = nexttile(t);
h = errorBarPlot(nanmean(behData.initAcc,2), 'type','bar');
ylim([.5 1]);
ylabel('mean acc');


% hist RT
ax = nexttile(t);
% hist([resp.time; resp.timeConf]', 0:.1:2);
ksdensities(col(behData.RT), 0:.1:2);
xlabel('RT');
ylabel('freq');
xlim([0 2]);

% score
ax = nexttile(t);
errorBarPlot(CountUnique(behData.score,2), 'type','bar');
set(gca,'XTick', 1:length(CountUnique(behData.score,2)), 'XTickLabel', unique(behData.score(~isnan(behData.score)))');
xlabel('score');
ylabel('freq');


% resp hand
ax = nexttile(t);
h = errorBarPlot(groupMeans(behData.initAcc, 2, behData.respLR), 'type','bar');
xticks(1:2);
xticklabels(labels.hand);
ylim([.5 1]);
ylabel('mean acc');
% legend(labels.hand,'Location','North');


% correct freq
ax = nexttile(t);
h = errorBarPlot(groupMeans(behData.initAcc,2, behData.corrFreq), 'type','bar');
xticks(1:2);
xticklabels(labels.freq);
ylim([.5 1]);
ylabel('mean acc');
 legend(labels.freq,'Location','North');


% amount of unique rt
ax = nexttile(t);
a = round(behData.RT*100);
a2 = unique (a);
a2 = a2(~isnan(a2));
 counts = histc (a(:), a2);
 bar (a2'/100, counts');
 xlabel ("RT");
ylabel ("amount of answers");

%% conf resps - errorbars across pps
% 
figure();
t = tiledlayout('flow');

% hist conf at both stages
ax = nexttile(t);
h = errorBarPlot(CountUnique(behData.certainty,2), 'type','bar');
h.FaceColor = 'flat';
h.CData = colours.certainty;
xticklabels(labels.certainty)
xlabel('certainty');
ylabel('freq');

% ksdensity RTs by factors - average over pps

figure();
t = tiledlayout('flow');

% rt by init acc
vars = {'RT'};

for i = 1:length(vars)
    ax = nexttile(t);
    x = 0:.1:2;
    f = permute(reshape(ksdensities(reshape(behDataByInitAcc.(vars{i}),nPP*2,[])', x),nPP,2,[]),[1,3,2]);
    if isfield(colours, 'initAcc'); set(gca,'ColorOrder',colours.initAcc,'nextplot','replacechildren'); end

    h = errorBarPlot(f, 'area',1,'xaxisvalues', x);

    xlabel(vars{i});
    ylabel('freq');
    xlim([0 2]);
    legend([h{:,1}], labels.initAcc,'Location','Best');
end

% rt by
vars2 = {'initResp','cue','certainty','confInCorr'};
for i = 1:length(vars2)
    ax = nexttile(t);
    y = groupMeans(behData.(vars{1}),2,behData.(vars2{i}),'dim');
    f = permute(reshape(ksdensities(reshape(y,nPP*size(y,2),[])', x),nPP,size(y,2),[]),[1,3,2]);
    if size(f,3)==6; set(gca,'ColorOrder',crameri('roma',6),'nextplot','replacechildren'); end
    if isfield(colours, vars2{i}); set(gca,'ColorOrder',colours.(vars2{i}),'nextplot','replacechildren'); end

    h = errorBarPlot(f, 'area',1,'xaxisvalues', x);


    %     ksdensities(reshape(permute(x,[1,3,2]),[],size(x,2)), 0:.1:2);
    xlabel('RT');
    ylabel('freq');
    xlim([0 2]);
    legend([h{:,1}], labels.(vars2{i}),'Location','Best');
end

%% making table with all 

%finding unique cues
unique_el = unique (behData.cue); unique_el = unique_el(1:3, :)';


%finding unique  correct/error
unique_acc = unique (behData.initAcc); unique_acc = unique_acc(1:2, :)';

%finding unique certainty
unique_cert = unique (behData.certainty); unique_cert = unique_cert(1:3, :)';


sorted_table = []; g =1; t_in_cell = [];

%loop for sorting evething based on cue (inv/n/valid), error/correct and
%certainty (maybe/prob/certain)

for i = 1:length (unique_el)
    for j = 1:length (unique_acc)
        for n = 1:length (unique_cert)
            c = sum (sum ( behData.certainty == unique_cert(n) & behData.cue == unique_el(i) & behData.initAcc == unique_acc(j)));
            t = behData.RT( behData.certainty == unique_cert(n) & behData.cue == unique_el(i) & behData.initAcc == unique_acc(j));
            t_m = mean (t);
            sd_t = std (t);
            t_in_cell{g} = t';
            pers_from_cue = length (t)/sum(sum (sum(behData.cue == unique_el(i))))*100;
            pers_from_cue_e = length (t)/sum(sum (sum(behData.cue == unique_el(i) & behData.initAcc == unique_acc(j))))*100;
            sorted_table = [sorted_table ; unique_el(i), unique_acc(j), unique_cert(n), c, t_m, sd_t, pers_from_cue, pers_from_cue_e, g];
            g = g+1;
        end
    end
end

t_in_cell =t_in_cell'; %timing in different cells

amount_of_each_answers = []; %Calculate the percentage of responses relative to the total number of occurrences
amount_of_each_answers = array2table ([sorted_table], 'VariableNames',{ 'cue', 'accuracy', 'certainty', 'amount of answers', 'mean timing of responses', 'sd time',...
    '% from this cue', '% from cue +corr/incorrect', '# of row in t_in_cell'}) %table with names for them

%making cells with labels
combined_labels = cell(length(unique_el) * length(unique_acc)*length(unique_cert), 1);
c = 1;

%making labels for graphs with RT
for i = 1:length(unique_el)
    for j = 1:length(unique_acc)
        for n = 1:length (unique_cert)
        combined_labels{c} = sprintf('%s %s %s', labels.cue{i}, labels.initAcc{j}, labels.certainty{n});
        c = c + 1;
    end
end
end 


%% Diff graphs based on privious table and cell array
x = 1:18; % for x-axes
%RT +cue + correct
% scatter with all RTs separated by cue, error/correct, 
figure;
for i = 1:18
    hold on
    scatter(i, t_in_cell{i}, 50, colors{i})  % Corrected 'flled' to 'filled'
end
hold off

xticks(1:length(amount_of_each_answers.("mean timing of responses")));
xticklabels(combined_labels);
xtickangle(60); % Rotating x-axis labels for better readability

% Labels and title
xlabel('Category');
ylabel('RT');
title('RT sorted by cue, correct-error and certainty');


%plot of sorted mean RT 
figure;
errorbar(x, amount_of_each_answers.("mean timing of responses"), amount_of_each_answers.("sd time"), 'o-', 'LineWidth', 2, 'MarkerSize', 8);
hold on
plot(x, polyval(polyfit(x, amount_of_each_answers.("mean timing of responses"), 1), x), '--');

xline (1:18, '--');
hold off

xticks(1:length(amount_of_each_answers.("mean timing of responses")));
xticklabels(combined_labels);
xtickangle(60); % Rotating x-axis labels for better readability

% Labels and title
xlabel('Category');
ylabel('RT');
title('RT sorted by cue, correct-error and certainty right against left hand');
legend ('mean RT');

%%%graphs of RT sorted by amount
% figure
% for i = 1 : length (t_in_cell)
% a = round(t_in_cell{i}*100);
% a2 = unique (a);
% a2 = a2(~isnan(a2));
%  counts = histc (a(:), a2);
% subplot (6, 3, i)
%  bar (a2'/100, counts', 'Facecolor', colors{i});
%  xlabel ( "RT for "+ combined_labels{i});
% ylabel ("amount of answers");
% 
% end

clear a a2

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
g = amount_of_each_answers.("# of row in t_in_cell") (amount_of_each_answers.("accuracy")==0);
k = amount_of_each_answers.("# of row in t_in_cell") (amount_of_each_answers.("accuracy")==1);

s_p =[];
for i = 1:9
    condition1 = t_in_cell{g(i)};
condition2 = t_in_cell{k(i)};
% Perform Friedman test
[p, ~, stats] = ranksum(condition1, condition2);
%disp(['p-value: ' num2str(p)]);
significant_points = p < 0.05;
s_p = [s_p, significant_points];
end
indices = find(s_p == 1);


%% errorbar 
y = 1:9;
figure
errorbar(y, amount_of_each_answers.("mean timing of responses")(amount_of_each_answers.("accuracy")==0), amount_of_each_answers.("sd time")(amount_of_each_answers.("accuracy")==0), 'o-', 'LineWidth', 2, 'MarkerSize', 6);
hold on
errorbar(y, amount_of_each_answers.("mean timing of responses")(amount_of_each_answers.("accuracy")==1), amount_of_each_answers.("sd time")(amount_of_each_answers.("accuracy")==1), 'o-', 'LineWidth', 2, 'MarkerSize', 6);

plot(y, polyval(polyfit(y, amount_of_each_answers.("mean timing of responses")(amount_of_each_answers.("accuracy")==0), 1), y), '--', 'LineWidth',1);
plot(y, polyval(polyfit(y, amount_of_each_answers.("mean timing of responses")(amount_of_each_answers.("accuracy")==1), 1), y), '--', 'LineWidth', 1);

plot(indices, zeros(size(indices)), 'r*');

set(gca,'ColorOrder',colours.initAcc,'nextplot','replacechildren')
xline (1:9, '--');
hold off

xticks(1: (length(amount_of_each_answers.("mean timing of responses")/2)));
xticklabels(combined_labels_2);
xtickangle(60); % Rotating x-axis labels for better readability

% Labels and title
xlabel('Category');
ylabel('RT');
title('RT sorted by cue and certainty');
legend ('RT for error',  'RT for correct', 'trend for error', 'trend for correct');


%% plot for mean RT for error and correct responce
%mean RT for errors
z  = 1:3;
figure
errorbar(z, amount_of_each_answers.("mean timing of responses")(amount_of_each_answers.("accuracy")==0 & amount_of_each_answers.("certainty")==1), amount_of_each_answers.("sd time")(amount_of_each_answers.("accuracy")==0 & amount_of_each_answers.("certainty")==1), 'o-', 'LineWidth', 2, 'MarkerSize', 8);
hold on
errorbar(z, amount_of_each_answers.("mean timing of responses")(amount_of_each_answers.("accuracy")==0 & amount_of_each_answers.("certainty")==2), amount_of_each_answers.("sd time")(amount_of_each_answers.("accuracy")==0 & amount_of_each_answers.("certainty")==2), 'o-', 'LineWidth', 2, 'MarkerSize', 8);
errorbar(z, amount_of_each_answers.("mean timing of responses")(amount_of_each_answers.("accuracy")==0 & amount_of_each_answers.("certainty")==3), amount_of_each_answers.("sd time")(amount_of_each_answers.("accuracy")==0 & amount_of_each_answers.("certainty")==3), 'o-', 'LineWidth', 2, 'MarkerSize', 8);
set(gca,'ColorOrder',colours.certainty,'nextplot','replacechildren')
xline (1:3, '--');
hold off

 xticks(1: (length(amount_of_each_answers.("mean timing of responses")/2)));
xticklabels(labels.cue);
xtickangle(60); % Rotating x-axis labels for better readability

% Labels and title
xlabel('Category');
ylabel('RT');
title('RT sorted by cue for error trials only');
legend ('RT for error + maybe', 'RT for error + probably', 'RT for error + certain');

%mean RT for correct
figure
errorbar(z, amount_of_each_answers.("mean timing of responses")(amount_of_each_answers.("accuracy")==1 & amount_of_each_answers.("certainty")==1), amount_of_each_answers.("sd time")(amount_of_each_answers.("accuracy")==1 & amount_of_each_answers.("certainty")==1), 'o-', 'LineWidth', 2, 'MarkerSize', 8);
hold on
errorbar(z, amount_of_each_answers.("mean timing of responses")(amount_of_each_answers.("accuracy")==1 & amount_of_each_answers.("certainty")==2), amount_of_each_answers.("sd time")(amount_of_each_answers.("accuracy")==1 & amount_of_each_answers.("certainty")==2), 'o-', 'LineWidth', 2, 'MarkerSize', 8);
errorbar(z, amount_of_each_answers.("mean timing of responses")(amount_of_each_answers.("accuracy")==1 & amount_of_each_answers.("certainty")==3), amount_of_each_answers.("sd time")(amount_of_each_answers.("accuracy")==1 & amount_of_each_answers.("certainty")==3), 'o-', 'LineWidth', 2, 'MarkerSize', 8);
set(gca,'ColorOrder',colours.certainty,'nextplot','replacechildren')
xline (1:3, '--');
hold off

 xticks(1: (length(amount_of_each_answers.("mean timing of responses")/2)));
xticklabels(labels.cue);
xtickangle(60); % Rotating x-axis labels for better readability

% Labels and title
xlabel('Category');
ylabel('RT');
title('RT sorted by cue for correct trials only');
legend ('RT for correct + maybe', 'RT for correct + probably', 'RT for correct + certain');


%% stats

% do pulses affect stuff?

% make a table, to do lme
alpha = .05;

behVarNames = {'initAcc','RT','certainty','confInCorr','chosenFreq','badTrial'};

nDV = length(behVarNames);
behDataNames = fieldnames(behData); % names of all behData

behNames = ['pp','sess','cue','respLR','corrFreq', behVarNames]; % not DVs
regTab = struct2table(structfun(@col, rmfield(behData, behDataNames(~ismember(behDataNames, behNames))), 'UniformOutput',0));
behNames =  regTab.Properties.VariableNames; % update


% make list of names with 'logistic' appended
logDVs = {'initAcc','badTrial'};
dvNames = behVarNames; % copy
logisticInds = ismember(behVarNames, logDVs);
dvNames(logisticInds) = strcat(dvNames(logisticInds), 'Logistic');

% make copy before nanzscoring
for j = find(logisticInds)
    regTab.(dvNames{j}) = regTab.(behNames{strcmp(behNames, behVarNames{j})});
end

% also get log RT
regTab.RTLog = log(regTab.RT);

behNames = regTab.Properties.VariableNames; % update

isLogistic = ismember(behNames, dvNames(logisticInds)); % which are [0 1] for logistic regressions

regTab(:,~isLogistic) = varfun(@nanzscore, regTab(:,~isLogistic));

glmeArgs = { {}, {'link','logit','distribution','binomial'} }; % for normal/logistic regression, if behVars are DV

% also get categorical cues, with neutral as reference
regTab.cueCateg = categorical(col(behData.cue),[0 -1 1], {'0','-1','1'});

ivNames = {'cue'}; % factors
nIV = length(ivNames);

%% main effects + RE

dvNames = {'initAccLogistic','RTLog','certainty','confInCorr','chosenFreq'};%,'badTrialLogistic'}; % not enough early/late yet
nDV = length(dvNames);

statNames = {'pValue','Estimate','Upper'};
stats = struct();
for ii = 1:length(statNames)
    stats.(statNames{ii}) = table();
end

fits = cell(nDV, nIV);
for j = 1:nIV
    for i = 1:nDV
        fits{i,j} = fitglme(regTab, sprintf('%s ~ 1 + %s + (1 | pp)', dvNames{i}, ivNames{j}), glmeArgs{logisticInds(i) +1}{:});
    end

    for ii = 1:length(statNames)
        stats.(statNames{ii}) = vertcat(stats.(statNames{ii}), StatsTableFromLMECells(fits(:,j), dvNames, statNames{ii}, 0));
    end
end

% get CI
stats.CI = array2table(table2array(stats.Upper) - table2array(stats.Estimate), 'VariableNames', stats.Estimate.Properties.VariableNames,'RowNames', stats.Estimate.Properties.RowNames);

figure();
c = get(gca,'ColorOrder');
set(gca,'ColorOrder', c(1:2,:), 'nextplot', 'replacechildren');
h = barwitherr(table2array(stats.CI)',table2array(stats.Estimate)');
xticks(1:nDV);
xticklabels(behVarNames);
ylabel(sprintf('Beta coefficient for the effect of %s', ivNames{1}));

hold on;
isSig = double(table2array(stats.pValue) <= alpha);
isSig(isSig==0) = NaN;
if size(isSig,1)>1
    plot(((1:nDV) +[-.1;.1])', isSig'*min(ylim), '*','LineWidth',2);
else
    plot(((1:nDV))', isSig'*min(ylim), '*','LineWidth',2);
end

if strcmp(ivNames, 'cueCateg')
    for i = 1:nDV
        isSig2(1,i) = double(anova(fits{i,j}).pValue(end) <= alpha);
    end
    
    isSig2(isSig2==0) = NaN;
    hold on;
    plot((1:nDV), isSig2'*min(ylim), 'k*','LineWidth',2);
end

title(sprintf('effect of %s', ivNames{1}));



%% making matrix with different conditions + corr for trails next to each other
 c = 1;
 all_con = [];
all_con = zeros(size (behData.block));

for i = 1:length (unique_el)
    for j = 1:length (unique_acc)
        for n = 1:length (unique_cert)
            a = find ( behData.certainty == unique_cert(n) & behData.cue == unique_el(i) & behData.initAcc == unique_acc(j));
            all_con(a) = c;
            c = c+1;
        end
    end
end

% matrix for cue, error and cert
many_cond = all_con';
many_cond = reshape (many_cond, [], 1);

all_con_1 = all_con';
m_RT = behData.RT';


%matrix for blocks
blocks = behData.pp';
blocks = reshape (blocks, [], 1);

%matrix for error
e_or_c = behData.initAcc';
e_or_c = reshape (e_or_c, [], 1);

%matrix for cue
cue_m = behData.cue';
cue_m = reshape (cue_m, [], 1);

%matrix for certainty
cert = behData.certainty';
cert = reshape (cert, [], 1);


%matrix for cue + error %%% 1 - inv cue + error, 2 - - inv cue + corr;
                        %%% 3 - n cue + error, 4 - - inv n + corr;
                        %%% 5 - valid cue + error, 6 - - valid cue + corr,
                        
c=1;
cue_and_e = [];
cue_and_e = zeros(size (behData.block));

for i = 1:3
    for j = 1:2
            a = find (  behData.cue == unique_el(i) & behData.initAcc == unique_acc(j));
            cue_and_e(a) = c;
            c = c+1;
  
    end
end

cue_and_e = cue_and_e';
cue_and_e= reshape (cue_and_e, [], 1);


%matrix for cue and cert %%% 1 - inv cue + maybe, 2 - inv cue + prob; 3 -inv cue + cert
                        %%% 4 - inv cue + maybe, 5 - inv cue + prob; 6 -inv cue + cert
                        %%% 7 - inv cue + maybe, 8 - inv cue + prob; 9 -inv cue + cert
 c = 1;
 cue_and_cert = [];
 cue_and_cert = zeros(size (behData.block));

for i = 1: 3
    for j = 1:3
            a = find (  behData.cue == unique_el(i) & behData.certainty == unique_cert(j));
             cue_and_cert(a) = c;
            c = c+1;
  
    end
end

cue_and_cert = cue_and_cert';
cue_and_cert= reshape (cue_and_cert, [], 1);


p = extractBefore(wh, '_');
pp = repmat({p}, length (cue_and_cert), 1);
%% saving file
fold = 'D:\cue_task\analysis\Data\Combine';
save (fullfile(fold, [wh '_sorted_t']),...
        'many_cond', 'cert',  'cue_and_cert', 'amount_of_each_answers','combined_labels', ...
        't_in_cell', 'behData',  'cue_and_e', 'e_or_c', 'cue_m', 'blocks');
