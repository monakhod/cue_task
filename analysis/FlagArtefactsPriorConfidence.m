% function FlagArtefactsPriorConfidence
% should flag from baseline until just after init resp?
% could try doing separate flaggings for cue-locked, ev-locked, stim-locked
% epochs

clear; close all; clc;
outFolder = 'D:\cue_task\analysis\Data\Saves';
loadUp = 0; % load the combined veog, isArt, etc that is already saved
useLong = 1; % flag to confRT + 250ms

if useLong 
    loadName = 'FlagArtefactsPriorConfidenceData_Long.mat'; % for initial loading of data
    saveName = 'FlagArtefactsPriorConfidence_Long.mat'; % for saving final file
    saveName1 = 'FlagArtefactsPriorConfidence_Long_';

% Concatenate strings

else
    loadName = 'FlagArtefactsPriorConfidenceData.mat'; % for initial loading of data
    saveName = 'FlagArtefactsPriorConfidence.mat'; % for saving final file
end

load(fullfile(outFolder, 'ExtractEpochsPriorConfidence.mat'));
fileInfo.interpFolder = 'D:\cue_task\analysis\Data\Interp\'; % my channels
ext = '_intp';


% for 1000Hz eye-tracker
art.rtLimsSamp = [102 1843]; % in samples [100 1800] ms
art.rtLimsTime = art.rtLimsSamp / eeg.fs * 1000; % in ms relative to evidence onset

eeg.veogChans = 129; % 131 is upper and 132 is lower - BUT CHECK THAT THESE ARE ACTUALLY YOUR VEOG CHANNELS - DIFFERENT PHD STUDENTS USE DIFFERENT CONVENTIONS!
art.blinkTh = 200; % threshold for detecting a blink in VEOG upper-lower difference signal (peak2peak moving)
art.blinkWindow = 200; % window (ms) for blink detection
art.blinkStep = 100; % stepsize (ms)
art.blinkTh2 = 100; % for eeglab func

art.saccTh = 50; % thresh for step-change (saccade detection)
art.saccWindow = 100; % ms
art.saccStep = 30; %ms

ppd = 25.9514; % 1 vis deg
art.ppd = ppd;
art.minSaccSizePix = ppd * 2;
art.targDist = ppd * 4; % max distance eye can go from start (baseline)

art.moveTh = [-200 200]; % [min max] change for detecting movement or other artefacts
art.artifTh = 100; % this is for SCALP channels - it's high because it is in CSD units, which take bigger values
art.chans = 1:2; % Artifact rejection channels - which electrodes should we reject based on?
art.windowTime = [-1400 2000]; % in msec, the limits of the time window relative to the event in which you will check for an artifact
% will need to change this to run up until final cue/resp?
art.windowInds= isBetween(eeg.epochTimes, art.windowTime); % indices of the timepoints for the artifact check window whose limits (min and max) are defined above

art.goodThresh = 1200; % min number good trials

[isBlink, nArtPerTr, isBadRT, isMissing, isTrigBlink, isTrigSacc,...
    isBlink2, isBlink3, isMovement, isOutFix, isOutTarg, isSaccade,...
    RT, isArt, acc, initResp] = deal(NaN(fileInfo.maxTr, fileInfo.nPP));
[maxArt, chanSDs] = deal(NaN(eeg.nChans, fileInfo.maxTr, fileInfo.nPP));
[trigsInWindows, saccAmplsInWindows] = deal(cell(fileInfo.nPP,1));
% erpArt must be overwritten each pp, due to size constraints
% put we can keep veog and the triggers
[dist,veog] = deal(NaN(fileInfo.nPP, sum(isBetween(eeg.epochTimes, art.windowTime)), fileInfo.maxTr));
dist = dist + NaN*1i; % make imag nan also (default is zero)
nTrs = NaN(fileInfo.nPP,1);

if loadUp && exist(fullfile(outFolder, loadName),'file')
    data = load(fullfile(outFolder, loadName),'fileInfo','dist',...
        'art','veog','nTrs','maxArt','isMovement','trigsInWindows','saccAmplsInWindows',...
        'RT','chanSDs','acc','initResp');

    % check loaded file matches the current script
    if ~equals(art, data.art)
        warning('current settings do not match loaded settings');
        % either change settings, do not loadUp, or overwrite data.art with art
        keyboard;
    end
    if ~equals(fileInfo, data.fileInfo)
        warning('fileInfo does not match');
        % may just be new pps added
%         keyboard;
    end


    ppsToDo = find(~ismember(fileInfo.ppID, data.fileInfo.ppID))'; % find any missings
%     ppsToDo = 8;
%     ppsToDo = 1:fileInfo.nPP; % or do all

    % need to change this, as newest won't always be at the end - would
    % need to re-size everything
    [new,old] = ismember(fileInfo.ppID, data.fileInfo.ppID);
    
    nT = size(data.dist,2);
    
    dist(new,1:nT,:) = data.dist(old(new),:,:); % move into new correct order
    veog(new,1:nT,:) = data.veog(old(new),:,:); 
    nTrs(new) = data.nTrs(old(new)); 
    maxArt(:,:,new) = data.maxArt(:,:,old(new)); 
    isMovement(:,new) = data.isMovement(:,old(new));
    RT(:,new) = data.RT(:,old(new)); 
    trigsInWindows(new) = data.trigsInWindows(old(new)); 
    saccAmplsInWindows(new,:,:) = data.saccAmplsInWindows(old(new)); 
    chanSDs(:,:,new) = data.chanSDs(:,:,old(new)); 
    acc(:,new) = data.acc(:,old(new));
    initResp(:,new) = data.initResp(:,old(new));


%     struct2workspace(data, 0); % do not overwrite

else
    ppsToDo = 1:fileInfo.nPP; % do all
end


% load up each indiv file, store veog, calc isArt + isMovement
if ~isempty(ppsToDo)
    for iPP = ppsToDo
        
        disp([fileInfo.ppID{iPP} '...'])
              
        data = load(fullfile(fileInfo.behFolder, [fileInfo.ppID{iPP} '_beh']), ... % load beh
            'RS', 'trigs','sTimes', 'saccAmpls',...
            'initAcc','initResp','respLR');
        data1 = load(fullfile(fileInfo.interpFolder, [fileInfo.ppID{iPP} ext]), 'erp'); % load EEG        
        data.erp=data1.erp; clear data1; % move to data
         
%         nTr = find(sq(sum(~isnan(data.erp),[1 2]))>0,1,'last');
        nTr = size(data.erp,3);
        nTrs(iPP) = nTr;
   
        
        erpArt = NaN(eeg.nChansTot2, sum(art.windowInds), nTr);
        [trigsInWindows{iPP}, saccAmplsInWindows{iPP}] = deal(cell(nTr,1));
        %% art window is -200: (RT + confRT + 100)
        for iTr = 1:nTr
            if ~isnan(data.RS(iTr))
   
    %             respLims = [art.windowTime(1) RS(iTr)/eeg.fs*1000]; % this goes from evOn-200 to initRT
                
                if useLong % confRT  + 250
                    respLims = [art.windowTime(1) (data.RS(iTr) +256)/eeg.fs*1000]; % this goes from baseline:RT+250ms
                else % conf RT
                    respLims = [art.windowTime(1) (data.RS(iTr))/eeg.fs*1000]; % this goes from baseline to RT
                end
    
                
                respInds = isBetween(eeg.epochTimes, respLims);
                erpArt(:,1:sum(respInds),iTr) = data.erp(1:eeg.nChansTot2, respInds, iTr); % get resp epoch
    
    
                % also store the triggers within that
                trigsInWindows{iPP}{iTr} = data.trigs{iTr}(isBetween(data.sTimes{iTr}, respLims));
                saccAmplsInWindows{iPP}{iTr} = data.saccAmpls{iTr}(isBetween(data.sTimes{iTr}, respLims));
            end
        end
        nT = size(erpArt,2);
        
        %% calc per pp
    
        maxArt(:,1:nTr,iPP) = sq(max(abs(erpArt(1:eeg.nChans,:,1:nTr)),[],2)); %[chans tr pp]
        isMovement(1:nTr,iPP) = sq(any(myArtextval(erpArt(1:eeg.nChansTot,:,:), art.moveTh))); %[chan 1 tr]. pop_artextval - min/max in window

        % also store SD over time per trial
        chanSDs(:,1:nTr, iPP) = nanstd(erpArt(1:eeg.nChans,:,1:nTr),[],2);

        dist(iPP,1:nT,1:nTr) = complex(erpArt(end-2,:,1:nTr), erpArt(end-1,:,1:nTr)); % eye-channels as complex (dist=abs)
        veog(iPP,1:nT,1:nTr) = erpArt(eeg.veogChans,:,1:nTr); % this is VEOG channel
    
        % store RT in ms
        RT(1:nTr,iPP) = data.RS ./ eeg.fs*1000; % (ms)

        % store acc + conf
        acc(1:nTr,iPP) = data.initAcc;
        initResp(1:nTr,iPP) = data.respLR;
    end

    % trim
    emptyT = all(isnan(veog) & isnan(dist), [1 3]);
    veog(:,emptyT,:) = [];
    dist(:,emptyT,:) = [];

    %% save that
    save(fullfile(outFolder, loadName),'fileInfo','dist',...
        'art','veog','nTrs','maxArt','isMovement','trigsInWindows','saccAmplsInWindows',...
        'RT','nTrs','acc','initResp','chanSDs');
end


%% now calc per pp - on the loaded up data
for iPP = 1:fileInfo.nPP
    nTr = nTrs(iPP);

    %% by trigger flags
    isTrigSacc(1:nTr,iPP) = cellfun(@(x) any(x >= art.minSaccSizePix), saccAmplsInWindows{iPP}); % start or end of blink
%     isTrigBlink(1:nTr,iPP) = cellfun(@(x) any(x==53 | x==56 | x==57), trigsInWindows{iPP}); % start or end of blink

    % saccade detection is the size of the saccade only
    % it does not account for any drifts that occur before it
    % so maybe apply some thresholding to see if eyes go outside of fixation at any point?
%     isOutFix(1:nTr,iPP) = any(abs(dist(iPP,:,1:nTr)) > art.minSaccSizePix,2);
    isOutTarg(1:nTr,iPP) = any(abs(dist(iPP,:,1:nTr)) > art.targDist); % outside targets

    %% try other approaches

    isArt(1:nTr,iPP) = sq(sum(maxArt(:,1:nTr,iPP) > art.artifTh)); % sum number of channels with abs max activity over thresh
    isBlink(1:nTr,iPP) = (max(veog(iPP,:,1:nTr),[],2)-min(veog(iPP,:,1:nTr),[],2)) > art.blinkTh;
    isBlink3(1:nTr,iPP) = max(abs(veog(iPP,:,1:nTr)),[],2) > art.blinkTh; % this will give a logical 1 if this threshold is exceeded on a certain trial, 0 otherwise
    isBlink2(1:nTr,iPP) = sq(myArtmwppth(veog(iPP,:,1:nTr), art.blinkTh2, art.blinkWindow, art.blinkStep, eeg.fs)); %[tr 1]
    isSaccade(1:nTr,iPP) = sq(myArtstep(veog(iPP,:,1:nTr), art.saccTh, art.saccWindow, art.saccStep, eeg.fs));        
    isMissing(1:nTr,iPP) = isnan(RT(1:nTr,iPP)) | sq(all(isnan(veog(iPP,:,1:nTr)),2));
    isBadRT(1:nTr,iPP) = ~isBetween(RT(1:nTr,iPP), art.rtLimsTime);
    
%     figure(1);clf;
%     bar(nansum(isArt(1:nTr,iPP),2)); xlabel('chan');ylabel('num arts');
%     drawnow;
%     %% plot some stuff
%     
%     clf;
%     
%     isOut = isBlink; %isBlink | sq(any(isArt,1)); %sq(any(isArt,1));
%     t = {'accept','reject'};
%     for i = 1:2
%         subplot(2,1,i);
%         plot(eeg.epochTimes, sq(veog(1,:,isOut(:,iPP) == i-1)));
%         ylim([-500 500]);
%         ylabel('\muV'); xlabel('time');
%         xline(-200); xline(2500); yline(0);
%         title(t{i});
%     end
    
end

%%  also need to check for average accuracy, and frequency of confidence responses

% get acc per sess + per block
accPerSess = nanmean(acc,1); % [1 nPP]
isLowAcc = accPerSess < .6;


freqs = struct();
% get freqs of confResps per sess + block
freqs.initResp = CountUnique(initResp,1); %[resp nPP]

% also get certainty
certainty = abs(round(initResp -3.5)); %1-6 -> [3 2 1 1 2 3]
freqs.certainty= CountUnique(certainty,1); %[cert nPP]

% normalise all, check if > .9
freqs = structfun(@(x) x ./ sum(x, 1), freqs, 'Uni',0);
freqChecks = structfun(@(x) double(any(x > .9, 1)), freqs, 'Uni',0);
isBadFreq = any(struct2array(freqChecks),3); %[1 nPP]

if any(isBadFreq) %|| any(isBadFreqBlock,'all')
    warning('bad frequency of confidences found, check them');
    keyboard;
end

%%
% isFlagged =  isTrigSacc==1 | isOutTarg==1 | isSaccade==1 | isBlink==1 | isMovement>0 ; % 1 = flagged for removal [tr pp]
isFlagged = isMissing==1 | isTrigSacc==1 | isOutTarg==1 | isBlink==1 | isMovement>0 | isBadRT==1 ; % 1 = flagged for removal [tr pp]
% isFlagged = isMissing==1 | isTrigSacc==1 | isOutTarg==1 | isBlink==1 | isArt>0 | isBadRT==1 | isBadConfRT==1; % 1 = flagged for removal [tr pp]

isGood = isFlagged==0;
isBad = isFlagged==1;


nGoodTrials = nTrs' - sum(isBad);
nBadTrials = sum(isBad);
ppIDs = cellfun(@(x) str2double(x(2:3)), fileInfo.ppID1);

if fileInfo.nPP > 1
    
    nGoodPerPP = groupMeans(nGoodTrials,2,ppIDs','dim')'; % [pp session]
    nBadPerPP = groupMeans(nBadTrials,2,ppIDs','dim')'; % [pp session]
    
    % disp([nGoodPerPP, sum(nGoodPerPP,2)]);
    % disp([nBadPerPP, sum(nBadPerPP,2)]);
    % disp([nGoodPerPP, sum(nGoodPerPP,2)]./ ([nGoodPerPP, sum(nGoodPerPP,2)] + [nBadPerPP sum(nBadPerPP,2)]));

else
    nGoodPerPP = nGoodTrials;
    nBadPerPP = sum(isBad);
end

disp(nansum(nGoodPerPP,2));
disp(round(nansum(nGoodPerPP,2) ./ nansum(nGoodPerPP + nBadPerPP, 2)*100,1));
% disp(nansum(nGoodTrials,'all'));
ppsToExclude = sum(nGoodPerPP,2) < art.goodThresh ; % [pp 1]
ppsToExcludeSess = col(repmat(ppsToExclude',4,1)); % per session

disp(nGoodPerPP);
disp(round(nGoodPerPP ./ (nGoodPerPP + nBadPerPP) * 100,1));

disp(nanmean(nGoodTrials ./ (nGoodTrials+nBadTrials))*100);
%% get each type separately

% get each one separately
allArtefacts = cat(3, isMissing, isTrigSacc, isOutTarg, isBlink, isMovement>0, isBadRT); %[tr pp art]
artTab = array2table(permute(nansum(allArtefacts),[2 3 1]),'VariableNames', {'Missing','Saccade','OutTarget','Blink','Art','RT1'},...
    'RowNames', fileInfo.ppID);
artTab.total = nBadTrials';
disp('number of artefactual trials for each artefact type:');
disp(artTab);


figure();
imagesc(table2array(artTab),[0 max(nTrs)]);
xticks(1:width(artTab));xticklabels(artTab.Properties.VariableNames);
yticks(1:fileInfo.nPP); yticklabels(fileInfo.ppID);
colorbar;
title('# artefacts per type');
   
% can I get unique artefacts for each type?
uniqArtTab = artTab;
fn = uniqArtTab.Properties.VariableNames;
nArts = repmat(nansum(allArtefacts,3),1,1,size(allArtefacts,3)); %[tr pp art]
uniqArtTab = array2table(permute(sum(allArtefacts==1 & nArts==1),[2 3 1]), 'VariableNames', artTab.Properties.VariableNames(1:end-1),'RowNames',fileInfo.ppID);
disp('number of unique artefactual trials for each artefact type:');
disp(uniqArtTab);


% %% plot num good trials
% figure();
% hist(nansum(nGoodPerPP,2));
% xline(art.goodThresh);
% xlabel('# good trials');
% ylabel('# pps');

%% plot each pp now

figure();
t = tiledlayout('flow');
for iPP = 1:fileInfo.nPP
    nexttile(t);
    xT = eeg.epochTimes(find(art.windowInds,1) + (0:size(veog,2)-1)); % times
    h = errorBarPlot(permute(groupMeans(veog(iPP,:,:),3,isGood(:,iPP),'dim'),[2,1,3]),'area',1,'xaxisvalues',xT);
    yline(0,'--k');
    if iPP==1;legend([h{:,1}], {'reject','include'},'Location','Best'); end
    ylabel('VEOG ');
    xlabel('ms from evidence start');
%     ylim([-1 1]);
    title(fileInfo.ppID(iPP));
end

%% plot veog - grand mean by incl/excl

if fileInfo.nPP > 1
    figure();
    xT = eeg.epochTimes(find(art.windowInds,1) + (0:size(veog,2)-1)); % times
%     veogByIncl = groupMeans(abs(dist./ppd),3,repmat(permute(isGood,[2,3,1]),1,size(veog,2))); %[pp t excl/incl]
    veogByIncl = groupMeans(veog,3,repmat(permute(isGood,[2,3,1]),1,size(veog,2))); %[pp t excl/incl]
    h = errorBarPlot(veogByIncl,'area',1,'xaxisvalues',xT);
    yline(0,'--k');
    legend([h{:,1}], {'reject','include'},'Location','Best');
    ylabel('VEOG');
    xlabel('ms from evidence start');
    
    %% plot all inlcuded together, sep pp lines
    figure();
    % veogByIncl = groupMeans(abs(dist./ppd),3,repmat(permute(isGood,[2,3,1]),1,size(veog,2)),'dim'); %[pp t excl/incl tr]
    veogByIncl = groupMeans(veog,3,repmat(permute(isGood,[2,3,1]),1,size(veog,2)),'dim'); %[pp t excl/incl tr]
    h = errorBarPlot(permute(veogByIncl(:,:,2,:),[4,2,1,3]), 'area', 1, 'xaxisvalues', xT);
    % h = plot(xT, veogByIncl(:,:,2)');
    ylabel('VEOG ');
    xlabel('ms from evidence start');
    title('seperate pp + sessions');
    legend([h{:,1}], fileInfo.ppID, 'Location','Best');
    
    %% plot all inlcuded together, sep pp lines, averaged across sessions (error bars are all trials in all sessions)
    
    figure();
    veogByInclPP = permute(groupMeans(veogByIncl, 1, ppIDs, 'dim'),[4,1,2,3,5]); %[pp, t, in/excl, sess, tr]
    veogByInclPP = reshape(veogByInclPP,size(veogByInclPP,1), length(xT), 2, []); %[pp, t, in/excl, sess*tr]
    h = errorBarPlot(permute(veogByInclPP(:,:,2,:),[4 2 1 3]), 'area', 1, 'xaxisvalues', xT);
    ylabel('VEOG ');
    xlabel('ms from evidence start');
    title('seperate pps (combined across sessions)');
    legend([h{:,1}], unique(fileInfo.ppID1), 'Location','Best');
end


%% plot indiv included trials per pp

figure();
n = 20; % number to show
t = tiledlayout('flow');
for iPP = 1:fileInfo.nPP
    if ~any(isGood(:,iPP)); continue; end
    nexttile(t);
    h = plot(xT, sq(abs(dist(iPP,:,find(isGood(:,iPP) & sq(~all(isnan(dist(iPP,:,:)),2)) ,n))))./ppd);
%     h = plot(xT, sq((veog(iPP,:,find(isGood(:,iPP) & sq(~all(isnan(dist(iPP,:,:)),2)) ,n)))));
    ylabel('eye-pos (vis deg)');
    xlabel('ms from evidence start');
%     ylim([0 4]);
    title(fileInfo.ppID(iPP));
end
title(t, 'example included trials');

%% plot eyepos averaged per pp

figure();
if fileInfo.nPP>1
    distByIncl = groupMeans(abs(dist./ppd),3,repmat(permute(isGood,[2,3,1]),1,size(veog,2)),'dim'); %[pp t excl/incl tr]
else
    distByIncl = permute(groupMeans(abs(dist./ppd),3,repmat(permute(isGood,[2,3,1]),1,size(veog,2)),'dim'),[4,1,2,3]); %[pp t excl/incl tr]
end
subplot(2,1,1); 
h = errorBarPlot(permute(distByIncl(:,:,2,:),[4,2,1,3]), 'area', 1, 'xaxisvalues', xT);
ylabel('Eye-distance (vis deg)');
xlabel('ms from evidence start');
title('seperate pp + sessions');
legend([h{:,1}], fileInfo.ppID, 'Location','Best');


%% combine across sessions
% 
% distByInclPP = permute(groupMeans(distByIncl, 1, ppIDs, 'dim'),[4,1,2,3,5]); %[pp, t, in/excl, sess, tr]
% distByInclPP = reshape(distByInclPP,size(distByInclPP,1), length(xT), 2, []); %[pp, t, in/excl, sess*tr]
% subplot(2,1,2);
% h = errorBarPlot(permute(distByInclPP(:,:,2,:),[4 2 1 3]), 'area', 1, 'xaxisvalues', xT);
% ylabel('Eye-distance (vis deg)');
% xlabel('ms from evidence start');
% title('seperate pps (combined across sessions');
% legend([h{:,1}], unique(fileInfo.ppID1), 'Location','Best');


%% plot eyepost

figure();
n = 20; % number to show
t = tiledlayout('flow');
for iPP = 1:fileInfo.nPP
    if ~any(isGood(:,iPP)); continue; end
    nexttile(t);
    h = plot(sq((dist(iPP,:,find(isGood(:,iPP) & sq(~all(isnan(dist(iPP,:,:)),2)) ,n))))./ppd);
    xlabel('x-coords eye-pos (vis deg)');
    ylabel('y-coords');
    hold on;
    DrawCircles([0 0], [.36 1 art.minSaccSizePix/ppd  4]'); % [fixRect, fix+arrows, 2degrees, targetRadius]
    axis([-4 4 -4 4]);
    axis('square');
    title(fileInfo.ppID(iPP));
end
title(t, 'example included trials');




%% need to have a quick look at overall channel activity/variance?
% i.e. are noisy trials still getting through?

% topoplot maxArt/chanSDs for accepted, per person, then overall

chanSDsGood = chanSDs; %[chan tr pp]
chanSDsGood(repmat(permute(isFlagged==1,[3 1 2]),eeg.nChans,1,1)) = NaN; % remove flagged
chanSDsGood = sq(nanmean(chanSDsGood,2)); %[chan pp]

mapLims = minMax(chanSDsGood,'all')';

% and mean
figure();
topoplot(nanmean(chanSDsGood,2),eeg.chanlocs,'electrodes','off','mapLimits',minMax(nanmean(chanSDsGood,2)));%,'colormap',cmap);
colorbar;
title('grand mean SD');

if fileInfo.nPP > 1
    figure();
    t = tiledlayout('flow'); ax=[];
    for i = 1:fileInfo.nPP
        ax(i) = nexttile(t);
        topoplot(chanSDsGood(:,i),eeg.chanlocs,'electrodes','off','mapLimits',mapLims);%,'colormap',cmap);
        colorbar;
        title(fileInfo.ppID{i});
    end
end
%% are they missing at random?

load(fullfile(outFolder, 'BehDataLoadPriorConfidence.mat'), 'behDataSess','labels');
% regression with means rejected by condition (+ acc?)

factors = {'initAcc','cues','certainty','RTLog','corrLR','day','block','trial'};
% for i = 1:length(factors)
%     try
%         means = groupMeans(isGood', 2, behData.(factors{i})); %[pp cond (interrupt/cont)]
%         anovas.(factors{i}) = rmanova(means, {'pp',factors{i}},'categorical',2);
%         disp(factors{i});
%         disp(anovas.(factors{i}));
%     end
% end

% make table - session data
behDataSess.RTLog = log(behDataSess.RT);
behTab = struct2table(structfun(@(x) nanzscore(col(x)), keepFields(behDataSess, [factors, 'pp1']),'Uni',0));


% regression
behTab.isGood = col(isGood);
formula = 'isGood ~ 1 + %s + (1 | pp1)';
fits = cell(1,length(factors));
for i = 1:length(factors)
    try
    fits{i} = fitglme(behTab, sprintf(formula, factors{i}), 'link','logit','distribution','binomial');
    end
end

stats = StatsTableFromLMECells(fits, factors);
disp(stats);
% seems that only accuracy differs, but not conf or pulses

% do a stepwise glm?
v = behTab.Properties.VariableNames;
pred = v(~ismember(v, {'isGood','pp','pp1','isBadTr','RTLog'})); % remove these
stepwise = stepwiseglm(behTab, 'isGood ~ 1 + initAcc + cues + certainty', 'ResponseVar', 'isGood', ...
            'PredictorVars', pred,...
            'link','logit','distribution','binomial');
% starting from constant, seems initAcc and confRTLog are sig pred
% same if some terms are included in the base model
disp(stepwise);

%% plot exclusions vs block/trial

figure();
t = tiledlayout('flow');

ax(1) = nexttile(t); % mean excl per block in session
errorBarPlot(groupMeans(isBad', 2, behDataSess.block),'type','line');
xlabel('block number in session');
ylabel('p(excluded)');

% cond plot exclusions per trial in session
ax(2) = nexttile(t); % mean excl per block in session
conditionalPlot(behDataSess.trial', isBad);
xlabel('trial number in session');
ylabel('p(excluded)');

% do it per block & miniblock
ax(3) = nexttile(t); % mean excl per block in session
miniBlockCount = behDataSess.miniBlock + behDataSess.block*2 - 2;
errorBarPlot(groupMeans(isBad', 2, miniBlockCount),'type','line');
xlabel('mini-block number in session');
ylabel('p(excluded)');

%% does mean eye-pos differ by corrLR?
% could this affect beta stuff, e.g. if eye normally moves in direction of
% dots, will give a lateralised potential towards chosen option? 
% use real(dist) to get x-distance

fac = 'initAcc'; % can change this to 'cues', 'initAcc', 'respLR', etc
labels.respLR = labels.hand;

distGood = dist ./ ppd; % convert to degrees
distGood(permute(repmat(isFlagged==1,1,1,size(dist,2)),[2 3 1])) = NaN; % remove flagged
if fileInfo.nPP > 1
    distByCorrLR = groupMeans(distGood, 3, permute(repmat(behDataSess.(fac),1,1,size(dist,2)),[1 3 2]), 'dim'); %[pp t l/r tr]
else
    distByCorrLR = permute(groupMeans(distGood, 3, permute(repmat(behDataSess.(fac),1,1,size(dist,2)),[1 3 2]), 'dim'),[4,1,2,3]); %[pp t l/r tr]
end

figure();
subplot(2,2,1);
h = errorBarPlot(nanmean(real(distByCorrLR),4), 'area', 1, 'xaxisvalues', xT);
ylabel('x-pos (deg)');
xlabel('ms from evidence start');
legend([h{:,1}], labels.(fac), 'Location','Best');
title('x-pos');

subplot(2,2,2);
h = errorBarPlot(nanmean(imag(distByCorrLR),4), 'area', 1, 'xaxisvalues', xT);
ylabel('y-pos (deg)');
xlabel('ms from evidence start');
legend([h{:,1}], labels.(fac), 'Location','Best');
title('y-pos');

subplot(2,2,3);
h = errorBarPlot(nanmean(abs(distByCorrLR),4), 'area', 1, 'xaxisvalues', xT);
ylabel('abs(distance) (deg)');
xlabel('ms from evidence start');
legend([h{:,1}], labels.(fac), 'Location','Best');
title('abs(dist)')

subplot(2,2,4); % use ITPC for angle
% h = errorBarPlot(angle(wITPC(distByCorrLR,4)), 'area', 1, 'xaxisvalues', xT);
h = errorBarPlot(nanmean(angle(distByCorrLR),4), 'area', 1, 'xaxisvalues', xT);
ylabel('angle ');
xlabel('ms from evidence start');
legend([h{:,1}], labels.(fac), 'Location','Best');
title('angle');
%% plotting what John said

% n = 11; % smoothing window
% figure();
% t = tiledlayout('flow');
% 
% ax = nexttile(t);
%  y = (groupMeans(isFlagged,1,behDataSess.block,'dim'));
% errorBarPlot(reshape(permute(reshape(y,6,72,2),[1 3 2]),12,72),'area',1)
% xline(60); xline(40); xline(20);
% 
% ax = nexttile(t);
% y = (groupMeans(isTrigSacc,1,behDataSess.block,'dim'));
% errorBarPlot(reshape(permute(reshape(y,6,72,2),[1 3 2]),12,72),'area',1)
% xline(60); xline(40); xline(20);
% 
% ax = nexttile(t);
% y = (groupMeans(isOutTarg,1,behDataSess.block,'dim'));
% errorBarPlot(reshape(permute(reshape(y,6,72,2),[1 3 2]),12,72),'area',1)
% xline(60); xline(40); xline(20);

%% save
disp (length (find (isGood(1:864)==1)))
 save(fullfile(outFolder, saveName), ...
    'isFlagged','isBlink','isArt','isTrigSacc','isTrigBlink',...
    'ppsToExclude','isBadRT', 'art','isMissing','nGoodTrials','isGood','nTrs',...
    'isBad','isMovement','isBlink2','isSaccade','isBlink3','RT',...
    'dist','veog','nBadTrials','ppsToExcludeSess','isOutTarg',...
    'artTab','uniqArtTab');

 fold = 'D:\cue_task\analysis\Data\Combine';
 s_name = cell2mat (fileInfo.ppNamesSess(1));
 s_name = [s_name '_flagged'];
 RS = data.RS';
 save(fullfile(fold, s_name), ...
    'isFlagged','isBlink','isArt','isTrigSacc','isTrigBlink',...
    'isGood', 'maxArt', 'nGoodTrials', 'RS');
% %% better SD check
% 
% % load up each file, epoch with same limits, exclude bad trials
% % take SD across trials rather than over time
% % then take 
% 
% [meanErp, chanTSD] = deal(NaN(fileInfo.nPP, sum(isBetween(eeg.epochTimes, art.windowTime)), eeg.nChans)); %[pp T chan]
% for iPP = 1:fileInfo.nPP
%     
%     disp([fileInfo.ppID{iPP} '...'])
%     
%     load(fullfile(fileInfo.interpFolder, [fileInfo.ppID{iPP} ext]), 'erp', 'RS','confRS')
%     
%     nTr = nTrs(iPP);
% 
%     
%     erpArt = NaN(eeg.nChansTot2, sum(art.windowInds), nTr);
%     [trigsInWindows{iPP}, saccAmplsInWindows{iPP}] = deal(cell(nTr,1));
%     % art window is -200: (RT + confRT + 100)
%     for iTr = 1:nTr
%         if ~isnan(RS(iTr)) && isFlagged(iTr,iPP)==0 % only do good trials
% 
% %             respLims = [art.windowTime(1) RS(iTr)/eeg.fs*1000]; % this goes from evOn-200 to initRT
%             
%             if useLong % confRT  + 250
%                 respLims = [art.windowTime(1) (RS(iTr) + confRS(iTr) +256)/eeg.fs*1000]; % this goes from evOn:confCue+250ms
%             else % conf RT
%                 respLims = [art.windowTime(1) (RS(iTr) + confRS(iTr))/eeg.fs*1000]; % this goes from evOn-200ms to confRT
%             end
% 
%             
%             respInds = isBetween(eeg.epochTimes, respLims);
%             erpArt(:,1:sum(respInds),iTr) = erp(1:eeg.nChansTot2, respInds, iTr); % get resp epoch
% 
% 
%             % also store the triggers within that
%             trigsInWindows{iPP}{iTr} = trigs{iTr}(isBetween(sTimes{iTr}, respLims));
%             saccAmplsInWindows{iPP}{iTr} = saccAmpls{iTr}(isBetween(sTimes{iTr}, respLims));
%         end
%     end
%     nT = size(erpArt,2);
%     
%     % calc per pp
% 
%     % also store SD across good trials
%     chanTSD(iPP,1:nT,:) = nanstd(erpArt(1:eeg.nChans,:,1:nTr),[],3)';
%     meanErp(iPP,1:nT,:) = nanmean(erpArt(1:eeg.nChans,:,1:nTr),3)';
% end
% 
% % trim
% chanTSD(:,all(isnan(chanTSD),[1 3]),:) = [];
% meanErp(:,all(isnan(meanErp),[1 3]),:) = [];
% 
% %% plot
% 
% % mean over time
% chanTSDMean = sq(nanmean(chanTSD,2))';  %[chan pp]
% 
% mapLims = minMax(chanTSDMean,'all')';
% figure();
% t = tiledlayout('flow'); ax=[];
% for i = 1:fileInfo.nPP
%     ax(i) = nexttile(t);
%     topoplot(chanTSDMean(:,i),eeg.chanlocs,'electrodes','off','mapLimits',mapLims);%,'colormap',cmap);
%     colorbar;
%     title(fileInfo.ppID{i});
% end
% 
% % and mean
% figure();
% topoplot(nanmean(chanTSDMean,2),eeg.chanlocs,'electrodes','off','mapLimits',minMax(nanmean(chanTSDMean,2)));%,'colormap',cmap);
% colorbar;
% title('grand mean SD');
% 
% %% plot over time too?
% 
% figure();
% imagesc(sq(nanmean(chanTSD))'); colorbar;
% title('mean SD');
% 
% %% plot means
% 
% 
% % mean over time
% meanErpSD = sq(nanstd(meanErp,[],2))';  %[chan pp]
% 
% mapLims = minMax(meanErpSD,'all')';
% figure();
% t = tiledlayout('flow'); ax=[];
% for i = 1:fileInfo.nPP
%     ax(i) = nexttile(t);
%     topoplot(meanErpSD(:,i),eeg.chanlocs,'electrodes','off','mapLimits',mapLims);%,'colormap',cmap);
%     colorbar;
%     title(fileInfo.ppID{i});
% end
% 
% % and mean
% figure();
% topoplot(nanmean(meanErpSD,2),eeg.chanlocs,'electrodes','off','mapLimits',minMax(nanmean(meanErpSD,2)));%,'colormap',cmap);
% colorbar;
% title('grand mean SD');
% 
% %% over time
% figure;
% h=errorBarPlot(meanErp,'area',1,'xaxisvalues',xT);
% xlim([-200 2000]);
% 
% % imagesc(sq(nanstd(meanErp))'); colorbar;
% 
