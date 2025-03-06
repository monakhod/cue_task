%sorting erp
clear all
clc

ERPfolder = 'D:\cue_task\analysis\Data\CSD'; %getting erp grom here
ARTfolder = 'D:\cue_task\analysis\Data\Combine'; %getting artifacts from here


allsubj = {'P21_1'};


filenameF = [ allsubj{1} '_flagged.mat'];  % Concatenate properly
filenameE = [allsubj{1} '_intp_csd.mat'];  % Concatenate properly
file_t = [allsubj{1} '_sorted_t.mat'];

load(fullfile(ERPfolder, filenameE)); 
load(fullfile(ARTfolder, filenameF)); 
load(fullfile(ARTfolder, file_t)); 


%% EEG params

eeg.fs = 1024; % sample rate
eeg.nChans = 128;  % number of EEG channels
eeg.nExts = 1;     % number of external channels (EOG etc)
eeg.etChans = 3;   % eye-tracking chans to add in, incl pupil
eeg.nChansTot = eeg.nChans + eeg.nExts;
eeg.nChansTot2 = eeg.nChans + eeg.nExts + eeg.etChans; % include Eye tracking channels

% load in a structure 'chanlocs' containing the x,y,z locations of each of the 128 scalp channels in the cap.
chanInfo = load('chanlocsBioSemi128.mat'); % note this was made by calling>> readlocs('cap128.loc') , and, since the locations for the 8 externals (129:136) are all (1,0,0), getting rid of those by calling chanlocs = chanlocs(1:128)
eeg.chanlocs = chanInfo.chanlocs;
eeg.chanNames = {eeg.chanlocs.labels}';

% define the epoch to extract (around ev onset)
eeg.epochLimsMS = [-1600 2100];    % note the use of integer number of cycles of SSVEP (17Hz)
eeg.epochSamples = round(eeg.epochLimsMS(1)/1000*eeg.fs):round(eeg.epochLimsMS(2)/1000*eeg.fs); % hence get the list of sample points relative to a given event that we need to extract from the continuous EEG
eeg.epochTimes = eeg.epochSamples*1000/eeg.fs; % hence get the timebase for the epoch (which we can plot average ERPs against) in milliseconds

%define the baseline window (vs evidence onset, which is 1200ms after cues)
eeg.blWin = [-1400 -1200]; % 

eeg.nSamples = length(eeg.epochSamples); % num timepoints in epoch

eeg.importEyeFlags = 1; % store eye-tracker saccade/blink/fix flags as triggers? (pop_importeyetracker)
eeg.detectEyeFlags = 0; % pop_detecteyemovements (>1 deg, ignores blink periods so none of those triggers)

% uncomment this if locking to confidence response
%  resp locking
eeg.confRespLimsMS = [-1600 250]; % around conf resp
eeg.confRespSamples = round(eeg.confRespLimsMS(1)/1000*eeg.fs):round(eeg.confRespLimsMS(2)/1000*eeg.fs); % hence get the list of sample points relative to a given event that we need to extract from the continuous EEG
eeg.confRespTimes = eeg.confRespSamples*1000/eeg.fs; % hence get the timebase for the epoch (which we can plot average ERPs against) in milliseconds
art.windowTime = [-1400 2000]; % in msec, the limits of the time window relative to the event in which you will check for an artifact
% will need to change this to run up until final cue/resp?
art.windowInds= isBetween(eeg.epochTimes, art.windowTime);

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
maxTr = size(erp,3);

colors = {[0, .447, .741], %blue
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
[0, 0, 0] %17 black
[0.3, 0.3, 0.3] %18 grey
[0.7, 0.7, 0.7]}; %19 lighter grey

%% mapping stimulus locked

k = find(eeg.epochTimes==0); %finding 0 point that corresponds for stim on set 

if ~exist('topoplot','file'); eeglab nogui; end % turning on eeg lab

%plotting map for window 0 to 200 ms (stim locked erp)
figure;
a = nanmean(erp (:,k:k+200,:), 3);
a = nanmean(a, 2);
topoplot(a,eeg.chanlocs,'electrodes', 'numbers',...
        'emarker2', {[], '.','k',10,1},'mapLimits',[-30 30]);%[0 max(nanmean(SD,2))]);
    colorbar;
    title([ num2str(300), 'ms - stim locked']);

%plotting maps for stim locked erp from 
figure();

for i = 1 :20
a = nanmean(erp (:,k,:), 3);
subplot(5,4,i);
    topoplot(a,eeg.chanlocs,'electrodes', 'on',...
        'emarker2', {[], '.','k',10,1},'mapLimits',[-30 30]);%[0 max(nanmean(SD,2))]);
    colorbar;
    title([ num2str(k-1639), 'ms ']);
    k = k+100;
end

%% making all matrices to pulling erp out
ImportantChans = [3 4 5 19 32];
isBadVolt = sq(any(maxArt(ImportantChans,:,:) > 100,1)); % [tr pp]
isBadVoltAll = sq(any(maxArt(:,:,:) > 200,1)); 
isGood (isBadVolt == 1) = 0;
isGood (isBadVoltAll == 1) = 0;
isGood = isGood(1:864);

%cleaning ERP
erp (:,:,isGood==0) = [];

%error-correct m
e_or_c (isGood==0) = [];

%cert
cert (isGood==0) = [];
 
%cue
cue_m (isGood==0) = [];

% cue +error-correct
cue_and_e (isGood==0) = [];

%cue+cert
 cue_and_cert (isGood==0) = [];

%18 varients out
many_cond (isGood==0) = [];

blocks (isGood==0) = [];
%checking amount of each answer
for i = 1:18
    s = sum (many_cond == i);
     fprintf('Amount of answers for %s is %s (before rejecting was %s).\n', char(combined_labels(i)), num2str (s), num2str (amount_of_each_answers.("amount of answers")(i)));
end 

%RT
RT = behData.RT';
RT = reshape (RT, [], 1);RT = RT*1000;
RT = RT*1000/eeg.fs;
RT(isGood==0) = [];


%% finding responce lock

t_for_each_resp = [];

for i = 1 : length (RT)
    h = find (eeg.epochTimes==250);
    tt = [(1: h)+round(RT(i))];
    t_for_each_resp = [t_for_each_resp; tt];
end

responce_matrix = [];
result_size = size (erp(:,1:length(eeg.confRespTimes), :));
responce_matrix = zeros(result_size);

for j = 1:length (RT)
resp_lock_erp = erp (:,t_for_each_resp(j,:), j);
responce_matrix(:,:,j) = resp_lock_erp;
end

%% mapping for resp locked
k = find(eeg.confRespTimes==0);
if ~exist('topoplot','file'); eeglab nogui; end

figure;
a = nanmean(responce_matrix (:,k-150:k-50,:), 3);
a = nanmean(a,2);
topoplot(a,eeg.chanlocs,'electrodes', 'numbers',...
        'emarker2', {[], '.','k',10,1},'mapLimits',[-30 30]);%[0 max(nanmean(SD,2))]);
    colorbar;
    title([ '150-50 ms before response ']);

figure();
k=k-400;
g=-400;
for i = 1 :12
 
a = nanmean(responce_matrix(:,k,:), 3);
subplot(3,4,i);
    topoplot(a,eeg.chanlocs,'electrodes', 'on',...
        'emarker2', {[], '.','k',10,1},'mapLimits',[-30 30]);%[0 max(nanmean(SD,2))]);
    colorbar;
    title([ num2str(g), 'ms (resp)']);
    k = k+50;
    g= g +50;
end

% k = find(eeg.confRespTimes==0);
% if ~exist('topoplot','file'); eeglab nogui; end
% 
% figure;
% a = nanmean(responce_matrix (:,k-150:k-50,cue_m ==0), 3);
% a = nanmean(a,2);
% topoplot(a,eeg.chanlocs,'electrodes', 'numbers',...
%         'emarker2', {[], '.','k',10,1},'mapLimits',[-30 30]);%[0 max(nanmean(SD,2))]);
%     colorbar;
%     title([ 'resp only for neutrals']);

% 
% k = find(eeg.confRespTimes==0);
% if ~exist('topoplot','file'); eeglab nogui; end
% 
% for i = 1:6
% figure;
% a = nanmean(responce_matrix (:,k-150:k-50,blocks == i), 3);
% a = nanmean(a,2);
% topoplot(a,eeg.chanlocs,'electrodes', 'numbers',...
%         'emarker2', {[], '.','k',10,1},'mapLimits',[-30 30]);%[0 max(nanmean(SD,2))]);
%     colorbar;
%     title(['resp for  block ' num2str(i)]);
% end


%applying important chennels
 names_chans = eeg.chanNames(ImportantChans);
 names_chan = eeg.chanNames(ImportantChans);
 names_chan = repelem(names_chan, 2);


%% making resp lock erp for different cond

%smooth erp
cpps_resp = movmean(responce_matrix, 102, 2);
cpps_stim  = movmean(erp, 102, 2);


figure;
errorBarPlot(permute(cpps_resp(ImportantChans,:,:), [3 2 1]), 'area', 1, 'alpha', 0.2,'xaxisvalues', eeg.confRespTimes)
xline(0)
yline(0)
title('resp locked erp for important chanels');
legend (names_chan);


%ploting resp lock erp for error/corr
a = nanmean (cpps_resp (ImportantChans, :, e_or_c ==0), 1); 
b = nanmean (cpps_resp (ImportantChans, :, e_or_c ==1), 1); 
figure;
errorBarPlot(permute(a(:,:,:), [3 2 1]), 'area', 1, 'alpha', 0.2,'xaxisvalues', eeg.confRespTimes)
hold on
errorBarPlot(permute(b(:,:,:), [3 2 1]), 'area', 1, 'alpha', 0.2,'xaxisvalues', eeg.confRespTimes, 'color', [0, 0.5, 0])
xline(0)
yline(0)
title(['resp locked erp for ' ' ' allsubj{1}]);
legend ('error', '', 'correct', '', 'Location','eastoutside');



%ploting resp lock erp for diff ccert
 a = nanmean (cpps_resp (ImportantChans, :, cert == 1), 1); 
 b = nanmean (cpps_resp (ImportantChans, :, cert == 2), 1); 
 c = nanmean (cpps_resp (ImportantChans, :, cert == 3), 1); 
figure;
errorBarPlot(permute(a(:,:,:), [3 2 1]), 'area', 1,'alpha', 0.2, 'xaxisvalues', eeg.confRespTimes, 'color', colors{19})
hold on
errorBarPlot(permute(b(:,:,:), [3 2 1]), 'area', 1,'alpha', 0.2, 'xaxisvalues', eeg.confRespTimes, 'color', colors{18})
errorBarPlot(permute(c(:,:,:),[3 2 1]), 'area', 1, 'alpha', 0.2,'xaxisvalues', eeg.confRespTimes, 'color', colors{17})
xline(0)
yline(0)
title(['resp locked erp for ' ' ' allsubj{1}]);
legend ('maybe', '', 'probably', '', 'certainty', '', 'Location','eastoutside');



%ploting resp lock erp for diff cues
 a = nanmean (cpps_resp (ImportantChans, :, cue_m ==-1), 1); 
 b = nanmean (cpps_resp (ImportantChans, :, cue_m ==0), 1); 
 c = nanmean (cpps_resp (ImportantChans, :, cue_m ==1), 1); 
figure;
errorBarPlot(permute(a(:,:,:), [3 2 1]), 'area', 1, 'alpha', 0.2, 'xaxisvalues', eeg.confRespTimes, 'color', colors{2})
hold on
errorBarPlot(permute(b(:,:,:), [3 2 1]), 'area', 1,'alpha', 0.2, 'xaxisvalues',  eeg.confRespTimes, 'color', colors{1})
errorBarPlot(permute(c(:,:,:),[3 2 1]), 'area', 1,'alpha', 0.2, 'xaxisvalues',  eeg.confRespTimes, 'color', colors{5})
xline(0)
yline(0)
title(['resp locked erp for ' ' ' allsubj{1}]);
legend ('invalid', '', 'neutral', '', 'valid', '', 'Location','eastoutside');



%ploting resp locked erp for cue+cert


%matrix for cue and cert %%% 1 - inv cue + maybe, 2 - inv cue + prob; 3 -inv cue + cert
                         %%% 4 - inv cue + maybe, 5 - inv cue + prob; 6 -inv cue + cert
                         %%% 7 - inv cue + maybe, 8 - inv cue + prob; 9 -inv cue + cert

k = 1;
for i = 1:3
figure ();
        a = nanmean (cpps_resp (ImportantChans, :, cue_and_cert == k), 1);
        b = nanmean (cpps_resp (ImportantChans, :, cue_and_cert == k+1), 1); 
        c = nanmean (cpps_resp (ImportantChans, :, cue_and_cert == k+2), 1);
        
        errorBarPlot(permute(a(:,:,:), [3 2 1]), 'area', 1, 'alpha', 0.2, 'xaxisvalues', eeg.confRespTimes, 'color', colors{2})
        hold on
        errorBarPlot(permute(b(:,:,:), [3 2 1]), 'area', 1,'alpha', 0.2, 'xaxisvalues',  eeg.confRespTimes, 'color', colors{1})
        errorBarPlot(permute(c(:,:,:),[3 2 1]), 'area', 1,'alpha', 0.2, 'xaxisvalues',  eeg.confRespTimes, 'color', colors{5})

xline(0)
yline(0)
title([labels.cue{i} ' resp locked erp for ' ' ' allsubj{1}]);
legend ('maybe', '', 'probably', '', 'certain','', 'Location','best');
hold off
k = k+(3);
end




%ploting resp lock erp for diff cues+error/correct
labels.cueInitAcc = {' error invalid', '',' error neutral', '',' error valid', '', ...
    ' correct invalid', '',' correct neutral', '',' correct valid', ''};

s=1
for i = 0:1
    figure;
 a = nanmean (cpps_resp (ImportantChans, :, e_or_c == i & cue_m ==-1), 1); 
 b = nanmean (cpps_resp (ImportantChans, :, e_or_c == i & cue_m ==0), 1); 
 c = nanmean (cpps_resp (ImportantChans, :, e_or_c == i & cue_m ==1), 1); 

errorBarPlot(permute(a(:,:,:), [3 2 1]), 'area', 1, 'alpha', 0.2, 'xaxisvalues', eeg.confRespTimes, 'color', colors {2})
hold on
errorBarPlot(permute(b(:,:,:), [3 2 1]), 'area', 1, 'alpha', 0.2, 'xaxisvalues', eeg.confRespTimes, 'color', colors {1})
errorBarPlot(permute(c(:,:,:),[3 2 1]), 'area', 1, 'alpha', 0.2, 'xaxisvalues', eeg.confRespTimes, 'color', colors {5})
xline(0)
yline(0)
title(['resp locked erp for ' labels.initAcc{i+1} ' ' allsubj{1}]);
legend (labels.cueInitAcc (s:s+5), 'Location','eastoutside');
s=6
end





%ploting resp lock erp for came cond (valid correct cert against all other)
for i = 1:3:18
    a = nanmean (cpps_resp (ImportantChans, :, many_cond == i), 1);
    b = nanmean (cpps_resp (ImportantChans, :, many_cond == i+1), 1);
    c = nanmean (cpps_resp (ImportantChans, :, many_cond == i+2), 1);
figure;
errorBarPlot(permute(a(:,:,:), [3 2 1]), 'area', 1 , 'alpha', 0.2, 'xaxisvalues', eeg.confRespTimes)
hold on
errorBarPlot(permute(b(:,:,:), [3 2 1]), 'area', 1, 'alpha', 0.2, 'xaxisvalues', eeg.confRespTimes, 'color', [0, 0.5, 0])
errorBarPlot(permute(c(:,:,:),[3 2 1]), 'area', 1, 'alpha', 0.2, 'xaxisvalues', eeg.confRespTimes, 'color', [0.8000 ,   0.2000  ,  0.5000])
xline(0)
yline(0)
title(['resp locked erp for '  ' ' allsubj{1}]);
legend (combined_labels{i}, '', combined_labels{i+1}, '', combined_labels{i+2}, '', 'Location','best');
end


%ploting resp lock erp for diff cond (valid correct cert against all other)
colors_1 = {[0.925, 0.6625, 0.549] %1  (lighter Orange (1))
    [.85, .325, .098], %2 orange - 2
    [0.425, 0.1625, 0.049] %3  (darker Orange (2))
    [0.925, 0.6625, 0.549] %1  (lighter Orange (1))
    [.85, .325, .098], %2 orange - 2
    [0.425, 0.1625, 0.049] %3  (darker Orange (2))

    [0.5, 0.724, 0.8705] %4 (lighter blue (1))
[0, .447, .741], %5 blue
         [0, 0.2235, 0.3705] %12  (darker blue (1))
         [0.5, 0.724, 0.8705] %4 (lighter blue (1))
[0, .447, .741], %5 blue
         [0, 0.2235, 0.3705] %12  (darker blue (1))

[0.733, 0.837, 0.594] %15  (lighter green (5))
          [.466, .674, .188], %green
[0.233, 0.337, 0.094] %16  (darker green (5))
[0.733, 0.837, 0.594] %15  (lighter green (5))
          [.466, .674, .188], %green
[0.233, 0.337, 0.094] %16  (darker green (5))
[0, 0, 0] %17 black
[0.3, 0.3, 0.3] %18 grey
[0.7, 0.7, 0.7]}; %19 lighter grey

for i = 1:3:18
    a = nanmean (cpps_resp (ImportantChans, :, many_cond == i), 1);
    b = nanmean (cpps_resp (ImportantChans, :, many_cond == i+1), 1);
    c = nanmean (cpps_resp (ImportantChans, :, many_cond == i+2), 1);
figure;
errorBarPlot(permute(a(:,:,:), [3 2 1]), 'area', 1,'alpha', 0.2,  'xaxisvalues', eeg.confRespTimes, 'color', colors_1{i})
hold on
errorBarPlot(permute(b(:,:,:), [3 2 1]), 'area', 1, 'alpha', 0.2, 'xaxisvalues', eeg.confRespTimes, 'color', colors_1{i+1})
errorBarPlot(permute(c(:,:,:),[3 2 1]), 'area', 1, 'alpha', 0.2, 'xaxisvalues', eeg.confRespTimes, 'color', colors_1{i+2})

xline(0)
yline(0)
title(['resp locked erp for ' ' ' allsubj{1}]);
legend ('',combined_labels{i}, '', combined_labels{i+1}, '', combined_labels{i+2},   'Location','best');
end






%% stim locked erp


%ploting stim lock erp for error/corr
% a = nanmean(cpps_stim(ImportantChans, :, e_or_c == 0), 1); 
% b = nanmean(cpps_stim(ImportantChans, :, e_or_c == 1), 1); 
% figure;
% errorBarPlot(permute(a(:,:,:), [3 2 1]), 'area', 1, 'alpha', 0.2, 'xaxisvalues', eeg.epochTimes)
% hold on
% errorBarPlot(permute(b(:,:,:), [3 2 1]), 'area', 1, 'alpha', 0.2, 'xaxisvalues', eeg.epochTimes, 'color', [0, 0.5, 0])
% xline(0)
% yline(0)
% title(['Stimulus locked ERP for ' ' ' allsubj{1}]);
% legend('error', '', 'correct', '', 'Location', 'eastoutside');
% 
% %ploting stim lock erp for diff ccert
% a = nanmean(cpps_stim(ImportantChans, :, cert == 1), 1); 
% b = nanmean(cpps_stim(ImportantChans, :, cert == 2), 1); 
% c = nanmean(cpps_stim(ImportantChans, :, cert == 3), 1); 
% figure;
% errorBarPlot(permute(a(:,:,:), [3 2 1]), 'area', 1, 'alpha', 0.2, 'xaxisvalues', eeg.epochTimes, 'color', colors{2})
% hold on
% errorBarPlot(permute(b(:,:,:), [3 2 1]), 'area', 1, 'alpha', 0.2, 'xaxisvalues', eeg.epochTimes, 'color', colors{1})
% errorBarPlot(permute(c(:,:,:), [3 2 1]), 'area', 1, 'alpha', 0.2, 'xaxisvalues', eeg.epochTimes, 'color', colors{5})
% xline(0)
% yline(0)
% title(['Stimulus locked ERP for ' ' ' allsubj{1}]);
% legend('maybe', '', 'probably', '', 'certainty', '', 'Location', 'eastoutside');



