function ExtractEpochsPriorConfidence(folderNames)
% Adapted from SimonPipeline/ExtractERPs.m and ExtractEpochsDelayedConf.m, by John Grogan, 2021.

if ~exist('pop_biosig','file')
    eeglab; close all;
end
close all;

%% set up folders

% change these folders to where you have saved the data
% fileInfo.dataFolder = {'../Task/Data/'}; % 1-16
fileInfo.dataFolder = {'D:\cue_task\analysis\Data\'}; % 1-16
fileInfo.rawFolder = 'D:\cue_task\analysis\Data\'; % where EEG saved
fileInfo.behFolder = 'D:\cue_task\analysis\Data\'; % where beh data saved
fileInfo.outFolder = 'D:\cue_task\analysis\Data\Saves'; % this is where the new files will be saved in later scripts

%% find files

if ~exist('folderNames','var')
%     folderNames = {
%                 'D02/*/';
%                 };
    p = 1;
    folderNames = cellfun(@(x) sprintf('P%02d', x), num2cell(p)', 'Uni',0);

end
% folders = dir(fileInfo.dataFolder{1});
% folderNames = {folders.name}';
% folderNames = folderNames(cellRegexpi(folderNames,'D2')>0); % remove some folders

f = [];
for i = 1:length(folderNames)
    f1 =  dir(fullfile(fileInfo.dataFolder{1}, folderNames{i}));
    [f1.thisFolder] = deal(folderNames{i}); % put this foldername in
    f = cat(1, f, f1);
end
fNames = {f.name}';

toKeep = cellRegexpi(fNames, '_\d+.bdf') > 0; % bdfs for experimental blocks

fNames(~toKeep) = []; % remove
f(~toKeep) = [];

% can exclude blocks here
toRemove = cellRegexpi(fNames, 'D15_2_2')>0;
fNames(toRemove) = [];
f(toRemove) = [];

if isempty(f)
    error('no files found');
end

ppNames = cellfun(@(x) x(1:regexp(x,'_')-1), fNames, 'UniformOutput', 0); % get ppNames
sesNums = cellfun(@(x) str2double(x(regexp(x,'P\d\d_')+4)), fNames, 'UniformOutput', 0); % get session number
sesNums = cat(1, sesNums{:});
blockNums = cellfun(@(x) str2double(x(regexp(x,'\.')-1)), fNames, 'UniformOutput',0); % block numbers
blockNums = cat(1, blockNums{:}); % remove from cells

ppNamesSess = cellfun(@(x) x(1:regexp(x,'_\d\.')-1), fNames, 'UniformOutput', 0); % include sesNum

% % uncomment this if combining across all ses
% if length(unique(sesNums(~isnan(sesNums))))>1 % adjust block numbers by session number
%     blockNums = blockNums + (sesNums-1)*6; % max 6 per session
% end

% days = arrayfun(@(x) str2double(x.folder(end)), f); % get day folder 
% 
% % append day2 blockNums - if laoding from Day1 + 2 folders
% blockNums(days==2) = blockNums(days==2) + 8;

fileInfo.files = f; % store
fileInfo.fNames = fNames;
fileInfo.ppNames = ppNames;
fileInfo.blockNums = blockNums;
fileInfo.sesNums = sesNums;
fileInfo.ppNamesSess= ppNamesSess;

fileInfo.ppID = unique(fileInfo.ppNamesSess,'stable'); % 1 per person_sess
% fileInfo.ppID = unique(fileInfo.ppNames,'stable'); % 1 per person_sess
fileInfo.ppID1 = cellfun(@(x) x(1:3), fileInfo.ppID, 'Uni',0); % 1 per pp
% fileInfo.ppID1n = cellfun(@(x) str2double(x(2:3)), fileInfo.ppID1); 
% fileInfo.ppID1n = col(repmat(1:12,4,1));
[~,~,fileInfo.ppID1n] = unique(fileInfo.ppID1);
fileInfo.nPP = length(fileInfo.ppID);

% maximums?
fileInfo.nBlocks = max(blockNums);
fileInfo.nTr = 240; % or 192?
fileInfo.maxTr = fileInfo.nBlocks * fileInfo.nTr;
% fileInfo.maxTr = 60;%max(sesNums); % if this is smaller than the number of trials a person has, then RS etc get filled with zeros on NaN trials

fileInfo.nTrs = repmat(fileInfo.nTr,1,fileInfo.nPP); % assume 720 per day for now
%% set triggers

triggers.exptStart = 36;
triggers.trialStart = 1;% cue appears
triggers.cueOff = [2 3 4]; % neutral, valid, invalid cue turns off, intro fade in stats
triggers.stimOnset = [8 9]; % baseline evidence starts, Left stim is at 20Hz or 25Hz
triggers.evOnset = [12 13]; % contrast difference starts, correct Answer is Left or Right;
triggers.resp = [26 27 28 29 30 31]; % [S D F J K L] conf+direction response
triggers.stimOffset = 15; % stim off
triggers.fb = 16;

triggers.blockBreak = [34 35]; % block break start/end
triggers.exptEnd = 40;

triggers.is20 = 8; % 20hz trigger
triggers.isLeft = 12; % [12 13]. whether corr ans is left

triggers.eyes = {'saccade$', 51; % must be end of string, to avoid clash with saccadeStart & saccadeEnd
                 'fixation', 52;
                 '_blink', 53;}; % L_blink (or R_blink) from eye-tracker, not from snipsaccades

% these are from the snipSaccades.m detection
triggers.eyes2 = {'saccadeStart', 54;
                  'saccadeEnd',   55;
                  'blinkStart',   56;
                  'blinkEnd',     57;  }; % will regexp the names, to deal with l/r

triggers.relevant = [1 2 3 4 5 8 9 12 13 15 16 26 27 28 29 30 31 34 35 36 40];

%% EEG params

eeg.reref = 0; % use average reference
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
% % conf resp locking
% eeg.confRespLimsMS = [-1500 500]; % around conf resp
% eeg.confRespSamples = round(eeg.confRespLimsMS(1)/1000*eeg.fs):round(eeg.confRespLimsMS(2)/1000*eeg.fs); % hence get the list of sample points relative to a given event that we need to extract from the continuous EEG
% eeg.confRespTimes = eeg.confRespSamples*1000/eeg.fs; % hence get the timebase for the epoch (which we can plot average ERPs against) in milliseconds

%% init filters

filters.doDetrend = 1;     % detrend the data (1) or not (0)?
filters.doLPF = 1;    % 1 = low-pass filter the data, 0=don't.
filters.doHPF = 0;    % high-pass filter the data? Usually we don't do this unless the drift in the data is bad and not linear (so detrending alone won't get rid of it)

filters.loCutOff = 40;       % Low Pass Filter cutoff in Hz
filters.hiCutOff = 0;    % Cutoff frequency for the HIGHPASS filter (which is a 'low cut off' of the overall banpass we're effectively doing). Typically 0.05 Hz and very rarely above 0.1 Hz because it can result in distortions of the ERPs - see Steve Luck stuff online

%% save this info

save(fullfile(fileInfo.outFolder, 'ExtractEpochsPriorConfidence.mat'),'fileInfo','filters','eeg','triggers');

%% extract and save each pp - looping through all available blocks

MEs = {};
for iPP = 1:fileInfo.nPP
%     try
        nTr(iPP) = extraction(iPP, fileInfo, eeg, filters, triggers);
%     catch ME
%         MEs = [MEs; ME];
%         disp(ME);
%     end

end

% %% get numbers of trials for everyone
% for iPP = 1:fileInfo.nPP
% % use this just to load and get nTr per each file
%     if exist(fullfile(fileInfo.rawFolder, [fileInfo.ppID{iPP} '_raw.mat']), 'file')
%         a = load( fullfile(fileInfo.rawFolder, [fileInfo.ppID{iPP} '_raw.mat']), 'RS');
%         nTr(iPP) = length(a.RS);
%     end
% end
% 
% 
% fileInfo.maxTr = max(nTr);
% fileInfo.nTrs = nTr;
% save(fullfile(fileInfo.outFolder, 'ExtractEpochsPriorConfidence.mat'),'fileInfo','-append');

end

function [iTr, epochs] = extraction(iPP, fileInfo, eeg, filters, triggers)
% epochs = extraction(iPP, fileInfo, eeg, filters, triggers)
% Loops through all blocks for 1 person, loads EEG, trims to relevant
% triggers, applies filters, extracts epochs, baselines, and gets some
% behavioural response data from epochs too. Saves outputs
iTr = NaN;
if exist(fullfile(fileInfo.rawFolder, [fileInfo.ppID{iPP} '_raw.mat']))
    fprintf('%s file already processed, skipping\n', [fileInfo.ppID{iPP} '_raw.mat']); 
    % return; 
end

%% extract

% Now loop through the subjects and extract the single trial ERPs and the other important informaiton about every trial:
    
    % should probably extend these per block?
    [corrLR, respLR, RS, blockNum, initResp, cues, certainty] = deal(NaN(1, fileInfo.maxTr));
    erp = NaN(eeg.nChansTot2, eeg.nSamples, fileInfo.maxTr); % incl eye-tracking chans
    blAmp = erp(:,1,:);
    [epochTrigs, epochSTimes,saccAmpls] = deal(num2cell(NaN(1, fileInfo.maxTr)));
    dataMat = NaN(fileInfo.maxTr, 18);
    nTrPerBlock = zeros(1,fileInfo.nBlocks);

    % now list the eeg files to be used for each subject:
    thisPP = find(cellRegexpi(fileInfo.fNames, fileInfo.ppID(iPP)) == 1);
    thisBlocks = fileInfo.blockNums(thisPP);
    
    fileNames = {};   
    for iB = 1:length(thisPP)%length(thisPP)-5:length(thisPP)%
        close all; % close eye sycn fig
        
        try % see if the file for this block exists?
            
            dataFolder = fileInfo.files(thisPP(thisBlocks(iB))).folder;
            fileNames{iB} = fullfile(dataFolder, fileInfo.files(thisPP(thisBlocks(iB))).name);
            fprintf('\nBlock %d: %s\n', fileInfo.blockNums(thisPP(thisBlocks(iB))), fileInfo.files(thisPP(thisBlocks(iB))).name);
            EEG = pop_biosig(fileNames{iB}); % read in EEG - this is an EEGLAB function that outputs a structure 'EEG' with lots of fileds, the most important of which is EEG.data - the actual EEG data!

        catch ME
            fprintf('\nfile not found: %s\n', fileNames{iB});
            disp(ME);
            continue; % skip to next block
        end

        if EEG.srate ~= eeg.fs 
            EEG = pop_resample(EEG,eeg.fs);
        end

        %%
        % remove empty chans
        veog = -diff(EEG.data(eeg.nChans + (1:2),:,:), [], 1); % make VEOG, will append (upper minus lower)
        EEG.data = EEG.data(1:eeg.nChans,:,:); % keep only EEG
        EEG.data(end+1,:,:) = veog; % store veog
        EEG.chanlocs = EEG.chanlocs(1:eeg.nChansTot);
        EEG.chanlocs(end).labels = 'VEOG';

        EEG.nbchan = size(EEG.data,1);
        
        % detrend?
        if filters.doDetrend
            EEG.data = detrend(EEG.data')';
        end
        
        if filters.doHPF
            [B,A] = butter(3,filters.hiCutOff/(eeg.fs/2),'high'); % butter inputs are order, cutoff freq (in normalised scale where 1=half sample rate, and optional third argument telling it to be e.g. 'high'-pass)
            EEG.data = filtfilt(B,A,double(EEG.data)')';
        end
        
        % LP Filter
        if filters.doLPF
            EEG.data = eegfilt(EEG.data,eeg.fs,filters.hiCutOff,filters.loCutOff);
        end 
        

        %% synchronise eye-traces into here

        f = what(dataFolder);
        pidStart = regexpi(fileNames{iB}, fileInfo.ppID{iPP}); % find ppID start
        pidStart = pidStart(end); % last one
        pp = fileNames{iB}(pidStart:end-4);
%         if strcmp(pp(end-1), '_') % if no leading zero on block num
%             pp = [pp(1:end-1) '0' pp(end)];
%         end

        elFileName = pp;
        edfName = [elFileName '.edf'];
        isEdf = exist(fullfile(dataFolder, [elFileName '.edf']),'file'); % does it exist?
        if ~isEdf % if cannot find edf
            % try old way - may fail if practice edf is stored e.g. P03_1.edf
    
            % find edf file name
            f = dir(dataFolder);
            f = {f.name}';
            
            isMatch = cellRegexpi(f, sprintf('_\\d?%s.edf', fileNames{iB}(end-4)))>0;
            if isMatch
                edfName = f{isMatch};
            elseif regexpi(fileNames{iB}, 'D02_1_1.bdf')
                disp('D02_1_1 edf is missing');
            else
                warning('edf file not found - manually search for it');
                keyboard;
            end
    
            % parse if not done
            elFileName = edfName(1:end-4); % just stem of name
            
            isEdf = exist(fullfile(dataFolder, [elFileName '.edf']),'file'); % update
        end
        
        if isEdf
            disp(edfName);

            % make an asc
            if ~exist( fullfile(dataFolder, [elFileName '.asc']), 'file' )
              edf2ascMat(fullfile(dataFolder, [elFileName '.edf']));
              if ~exist( fullfile(dataFolder, [elFileName '.asc']), 'file' ) % if still not exist, because edf2asc won't run on one-Drive
                  edf2ascLoopPriorConfidence({fileInfo.files(thisPP(thisBlocks(iB))).thisFolder}); % will copy to local, convert, copy back
              end
            end
            
            
            if ~exist( fullfile(dataFolder, [elFileName '_el.mat']), 'file' )
              ET = parseeyelink( fullfile(dataFolder, [elFileName '.asc']), fullfile(dataFolder, [elFileName '_el.mat']),'TRIG');
            else
              ET = load( fullfile(dataFolder, [elFileName '_el.mat']) );
            end
    
            synchTrigs = [36 40];
            % some blocks are missing final trigger (endexpt), and instead
            % end on trigger 34, so use those
            % record the names of these blocks here
            blocksMissingEndTrigs = {['P14_1_1']};
    
            if any(cellfun(@any, regexp(fileNames{iB}, blocksMissingEndTrigs))) && (strcmp(EEG.event(end).type,'34') || EEG.event(end).type==34) 
                % use trig 34 as final synchronising trigger
                synchTrigs(2) = 34;

            elseif regexp(fileNames{iB}, 'D10_2_6')  %I'll leave this here to show how to deal with it
                % EEG recording stopped halfway through trial 111 (trig 15)
                % rename that here to trig40. I have manually changed it in
                % the asc file to trig40 (and trig40->trig41) too
                EEG.event(end).type = '40';
                EEG.event(end).edftype = [];

                ET.event(996:end,:) = []; % remove those additional triggers from ET
    
            elseif ~strcmp(EEG.event(end).type, '40')
                % the recording was probably stopped before the final 'trig40' was sent
                % so it probably ends on trigger 34. if so, add the block
                % name to 'blocksMissingEndTrigs' above and re-run.
                % If it is something else, you can change/remove triggers
                % from EEG to make them match, and/or rename triggers in asc
                
                keyboard;
%                 if strcmp(EEG.event(end).type, '34')
%                     
%                     % recording was stopped before 'end expt' trigger, so
%                     % just use final block break trigger
%                     synchTrigs = [36 34];
%                 end
    
            end
           
            % import that mat file into the EEG file - don't keep TIME or PUPIL. 
            % make input after labels 1 to store saccades as triggers
            EEG = pop_importeyetracker(EEG, fullfile(dataFolder, [elFileName '_el.mat']),synchTrigs ,[2 3 4] ,{'elX' 'elY' 'pupil'},eeg.importEyeFlags,1,0,1,10);
            % doesn't seem to store all blinks and saccades?
            drawnow;
    
            if max(abs(EEG.etc.eyetracker_syncquality(EEG.etc.eyetracker_syncquality(:,2)~=0,1))) > 1 || sum(EEG.etc.eyetracker_syncquality(:,2))<100
                disp(EEG.etc.eyetracker_syncquality(EEG.etc.eyetracker_syncquality(:,2)~=0,:));
                warning('please check synchronisation');
    %             keyboard;
                WaitSecs(2);
            end
            % maybe automate this?
                
            if eeg.detectEyeFlags % detects saccades over 1 degree
                
                if eeg.importEyeFlags % if already imported them, remove all the fix/saccades, but keep blinks
                    % will then detect saccades >1deg below
    
                    inds = cellRegexpi({EEG.event.type}', 'saccade|fixation')>0;
                    EEG.event(inds) = [];
    
                    % and from urevent?
                end
    
                EEG = pop_detecteyemovements(EEG, eeg.nChansTot + (1:2), [],6,20, [],1,0, 51, 2,1,1,1);
                drawnow;
            end
        
    
        
            %% snipsaccades to detect saccades/blinks?
    
            % inserts trigger codes for start/end of saccades/blinks into
            % EEG.event. if minSaccSizePix=0, does all saccades, otherwise set 
            % minimum saccade size (in pixels) to include
            
            minSaccSizePix = 0; % include all ; 26.4773; % only >1 vis deg
            [EEG, triggers.eyes2] = InsertSaccadesBlinksAsTriggers(EEG, elFileName, dataFolder, minSaccSizePix, triggers.exptStart); % get all saccades, will flag by ampltiude later

        else
            warning('edf file missing, so skipping all edf synchronisation');
            
            % need to manually insert blank channels to match missing ones?
            EEG.data(130:132,:) = NaN;
            EEG.nbchan = size(EEG.data,1);
            % add in channel names too
            elNames  = {'elX', 'elY','pupil'};
            for i = 1:3
                EEG.chanlocs(end+i).labels = elNames{i};
                EEG.chanlocs(end).types = 'EYE';
            end
            clear elNames;

            % need to also create 'sac_amplitude' field in EEG.event
            for i = 1:length(EEG.event)
                EEG.event(i).sac_amplitude = 0;
            end
        end


        %% load beh data

        f = what(dataFolder);
        pidStart = regexpi(fileNames{iB}, fileInfo.ppID{iPP}); % find ppID start
        pidStart = pidStart(end); % last one
        pp = fileNames{iB}(pidStart:end-4);
%         if strcmp(pp(end-1), '_') % if no leading zero on block num
%             pp = [pp(1:end-1) '0' pp(end)];
%         end
        try % TEST
            matName = f.mat{cellRegexpi(f.mat, sprintf('%s_Test', pp))>0};
        catch % ERROR
            matName = f.mat{cellRegexpi(f.mat, sprintf('%s-Error', pp))>0};
        end
%         matName = [fileNames{iB}(1:end-4) '.mat'];
        mat = load(fullfile(dataFolder, matName));


        % this is order of columns in matrix I will make
        dataMatNames = {'pp','day','block','miniBlock','deltaC',...
            'cues','trialLR','gratingFreq','targetFreq', 'stimShift', ...
            'initLR','initRT','initAcc','initResp','certainty',...
            'chosenFreq', 'reward','isBadTrial'};

        nTr = length(mat.par.trialLR); % not accurate, as are pre-filled
        nTrPerBlock(iB) = nTr; % append couter

        % some fields have smaller size, due to replay issue
        if any(structfun(@length, mat.resp) ~= nTr)
            fn = fieldnames(mat.resp);
            fInds = find(structfun(@length, mat.resp) ~= nTr); % ones mismatched
            for j = 1:length(fInds)
                nn = nTr - length(mat.resp.(fn{fInds(j)})); % fill with NaN
                if ~iscell(mat.resp.(fn{fInds(j)}))
                    mat.resp.(fn{fInds(j)}) = [mat.resp.(fn{fInds(j)}), NaN(size(mat.resp.(fn{fInds(j)}),1), nn)];
                else
                     mat.resp.(fn{fInds(j)}) = [mat.resp.(fn{fInds(j)}), cell(size(mat.resp.(fn{fInds(j)}),1), nn)];
                end

            end
        end

        % get them now, 1 row per trial (incl all trials in this block)
        origMat = [repmat(iPP,1,nTr); repmat(mat.par.day, 1,nTr);
            repmat(iB,1,nTr); mat.resp.block(1:nTr);
            repmat(mat.par.deltaC,1,nTr); 
            mat.par.cues(1:nTr); mat.par.trialLR(1:nTr);
            mat.par.gratingFreq(1:nTr); mat.par.targetFreq(1:nTr);
            mat.par.stimShift(1:nTr);
            mat.resp.LR(1:nTr); mat.resp.time(1:nTr)*1000; 
            mat.resp.perf(1,1:nTr); mat.resp.initResp(1:nTr); 
            mat.resp.certainty(1:nTr); 
            mat.resp.freq(1:nTr); mat.resp.reward(1:nTr); 
            mat.resp.badTrial(1:nTr);
            ]';

        %% find triggers

        % Fish out the event triggers and times
        nEvents = length(EEG.event); % total number of events
        
        if nEvents == 0
            warning('no trigger events found in this block (duration is %d), skipping', length(EEG.times)/eeg.fs);
            continue; 
        end % no triggers in this block?
        
        trigs = [EEG.event.edftype]; % get number triggers
        if length(trigs) < nEvents % if fewer triggers found than expected
            emptyInds = find(cellfun(@isempty, {EEG.event.edftype})); % find missing ones
            for i = 1:length(emptyInds) % convert them from the string in .type
                EEG.event(emptyInds(i)).edftype = str2double(EEG.event(emptyInds(i)).type);
            end

            trigs = [EEG.event.edftype]; % get them again
            if length(trigs) < nEvents % if there are still missing ones

                trigs = {EEG.event.type}; % this is sometimes 'condition X' or just 'X' (where X is a string number)
                if iscellstr(trigs)
                    % get triggers as numbers only
                    trigsNum = NaN(1,nEvents);
                    for i = 1:nEvents
                        trigsNum(i) = str2double(trigs{i}(regexp(trigs{i}, '\d')));
                    end
                    trigs = trigsNum; % replace
                elseif all(isnumeric([trigs{:}]))
                    trigs = [trigs{:}]; % 
                else
                    keyboard;
%                     error('trigs not converted to numbers');
                end
            end
        end

        sTimes = round([EEG.event.latency]);

        % so trigs and stimes are the exact same length; trigs has the trigger codes and stimes the SAMPLE POINTS in the
        % continuous EEG where those triggers happened.

        % some people had wrong trigger codes - they had 768 added for some reason?
        % check and fix
        if ~any(ismember(trigs, triggers.relevant))
            fprintf('\n no relevant trigger codes found, will try adjusting by 768 to fix known error...')
            keyboard;
            if all(ismember(trigs-768, [triggers.relevant 0]))
                trigs = trigs - 768; % minus
                trigs(trigs==0) = 1; % the 'start' should be 1 not zero
                if all(ismember(triggers.relevant, trigs))
                    fprintf(' fixed!\n');
                else
                    fprintf(' NOT FIXED - manually adjust\n');
                    beep;
                    keyboard;
%                     error('wrong trigs found');
                end
            else
                fprintf(' NOT FIXED - manually adjust\n');
                beep;
                keyboard;
%                 error('wrong trigs found')
            end
        end
        
        if ~any(trigs(1) == [36 1 0])
            trigs(1)=36; % it's usually because the recording started late, so set this to 'start of 
            keyboard;
%             error('wrong first trig');
        end

        % re-label eye-tracker trigger flags
        if eeg.importEyeFlags
            trigNames = {EEG.event.type}';

            for i = 1:size(triggers.eyes,1)
                inds = cellRegexpi(trigNames, triggers.eyes{i,1})>0; % find matches

                % check those triggers are presently zero
                if ~all(trigs(inds) ==0)
%                     error('some eye-tracking trigger labels are not zero');
                    keyboard;
                
                end

                % replace
                trigs(inds) = triggers.eyes{i,2};
            end
        end

        %% check all triggers match up

        if isEdf
            eegRelTrig = ismember(trigs, triggers.relevant);
            etRelTrig = ismember(ET.event(:,2), triggers.relevant);
            [freq, tr] = CountUnique(trigs(eegRelTrig)');
            [freq2, tr2] = CountUnique(ET.event(etRelTrig,2));
            if sum(eegRelTrig) ~= sum(etRelTrig) || ~equals(freq, freq2)
                warning('different number of relevant triggers in EEG and ET');
    
    
                if ~equals(tr, tr2)
                    warning('entire trigger missing from one');
                    if synchTrigs(2) == 34 && length(tr2)>length(tr) && tr2(~ismember(tr2,tr))==40
                        warning('was just due to stopping eeg recording a bit too early and missing trig40');
                    end
                elseif ~equals(freq, freq2)
                    warning('mismatch in triggers below (4th column):')
                    try; disp([tr, freq, freq2, freq~=freq2]); end
                    keyboard;
                end
                
                ind = find(diff(nancat(2, trigs(eegRelTrig)', ET.event(etRelTrig,2)),[],2),1); % first one to mismatch
                % find ind in triggers
                eegInd = find(eegRelTrig, ind); eegInd = eegInd(end); % get ind-th one
                etInd = find(etRelTrig, ind); etInd = etInd(end); % get ind-th one
    
                % display +/-3
                try
                disp(trigs(eegInd + (-3:3)));
                disp(ET.event(etInd + (-3:3),2)');
                end
%3 block
                if regexp(fileNames{iB}, 'P15_3_1') % can fix known issues like this here
                    if ind==801 && trigs(eegInd)==160 
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd) = []; % delete
                    % update these
                    trigs(eegInd) = [];
                    sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end

                    end

                    elseif regexp(fileNames{iB}, 'P20_1_4') % can fix known issues like this here
                    if ind==153 && trigs(eegInd-1)==157
                    % trig 13 missing, seems replaced by trg 45 (times
                    % match ET exactly, so overwriting it)
                    EEG.event(eegInd*1).type = []; % delete
                    % update these
                     EEG.event(eegInd-1).edftype = []; %
                    trigs(eegInd-1) = [];
                     sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end
                    elseif regexp(fileNames{iB}, 'P20_3_5') % can fix known issues like this here
                    if ind==70 && trigs(eegInd)==219
                    % trig 13 missing, seems replaced by trg 45 (times
                    % match ET exactly, so overwriting it)
                    EEG.event(eegInd).type = []; % delete
                    % update these
                     EEG.event(eegInd).edftype = []; %
                    trigs(eegInd) = [];
                     sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end

                     elseif regexp(fileNames{iB}, 'P20_3_6') % can fix known issues like this here
                    if ind==584 && trigs(eegInd)==223
                    % trig 13 missing, seems replaced by trg 45 (times
                    % match ET exactly, so overwriting it)
                    EEG.event(eegInd).type = []; % delete
                    % update these
                     EEG.event(eegInd).edftype = []; %
                    trigs(eegInd) = [];
                     sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end

 elseif regexp(fileNames{iB}, 'P22_2_5') % can fix known issues like this here
                    if ind==584 && trigs(eegInd)==223
                    % trig 13 missing, seems replaced by trg 45 (times
                    % match ET exactly, so overwriting it)
                    EEG.event(eegInd).type = []; % delete
                    % update these
                     EEG.event(eegInd).edftype = []; %
                    trigs(eegInd) = [];
                     sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end

                     elseif regexp(fileNames{iB}, 'P21_2_5') % can fix known issues like this here
                    if ind==569 && trigs(eegInd)==31
                    % trig 13 missing, seems replaced by trg 45 (times
                    % match ET exactly, so overwriting it)
                    EEG.event(eegInd).type = []; % delete
                    % update these
                     EEG.event(eegInd).edftype = []; %
                    trigs(eegInd) = [];
                     sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end


                     elseif regexp(fileNames{iB}, 'P15_3_2') % can fix known issues like this here
                    if ind==29 && trigs(eegInd)==191
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd).type = []; % delete
                    % update these
                     EEG.event(eegInd).edftype = []; %
                    trigs(eegInd) = [];
                     sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end

                     elseif regexp(fileNames{iB}, 'P15_3_3') % can fix known issues like this here
                    if ind==902 && trigs(eegInd)==40
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd).type = []; % delete
                    % update these
                     EEG.event(eegInd).edftype = []; %
                    trigs(eegInd) = [];
                     sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end

                     elseif regexp(fileNames{iB}, 'P15_3_4') % can fix known issues like this here
                    if ind==423 && trigs(eegInd)==35
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd).type = 3; % delete
                    % update these
                     EEG.event(eegInd).edftype = 3; %
                    trigs(eegInd) = 3;
                     %sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end


                     elseif regexp(fileNames{iB}, 'P15_3_5') % can fix known issues like this here
                    if ind==558 && trigs(eegInd)==33
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd).type = []; % delete
                    % update these
                     EEG.event(eegInd).edftype = []; %
                    trigs(eegInd) = [];
                     sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end
        elseif regexp(fileNames{iB}, 'P15_3_5') % can fix known issues like this here
                    if ind==605 && trigs(eegInd-1)==187
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd-1).type = []; % delete
                    % update these
                     EEG.event(eegInd-1).edftype = []; %
                    trigs(eegInd-1) = [];
                     sTimes(eegInd-1) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end
  %1 block
 elseif regexp(fileNames{iB}, 'P15_1_1') % can fix known issues like this here
                    if ind==358 && trigs(eegInd)==191
                    % trig 13 missing, seems replaced by trg 45 (times
                    % match ET exactly, so overwriting it)
                    EEG.event(eegInd).type = 3; % delete
                    % update these
                     EEG.event(eegInd).edftype = 3; %
                    trigs(eegInd) = 3;
                    % sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end
                   
 elseif regexp(fileNames{iB}, 'P15_1_1') % can fix known issues like this here
                    if ind==386 && trigs(eegInd-1)==186
                    % trig 13 missing, seems replaced by trg 45 (times
                    % match ET exactly, so overwriting it)
                    EEG.event(eegInd*1).type = []; % delete
                    % update these
                     EEG.event(eegInd-1).edftype = []; %
                    trigs(eegInd-1) = [];
                     sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end
 elseif regexp(fileNames{iB}, 'P15_1_2') % can fix known issues like this here
                    if ind==150 && trigs(eegInd)==33
                    % trig 13 missing, seems replaced by trg 45 (times
                    % match ET exactly, so overwriting it)
                    EEG.event(eegInd).type = []; % delete
                    % update these
                     EEG.event(eegInd).edftype = []; %
                    trigs(eegInd) = [];
                     sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end

                     elseif regexp(fileNames{iB}, 'P15_1_4') % can fix known issues like this here
                    if ind==863 && trigs(eegInd)==47
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd).type = []; % delete
                    % update these
                     EEG.event(eegInd).edftype = []; %
                    trigs(eegInd) = [];
                     sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end

                    %2 block

                         elseif regexp(fileNames{iB}, 'P15_2_5') % can fix known issues like this here
                    if ind==214 && trigs(eegInd)==32
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd).type = []; % delete
                    % update these
                     EEG.event(eegInd).edftype = []; %
                    trigs(eegInd) = [];
                     sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end

                    elseif regexp(fileNames{iB}, 'P15_2_1') % can fix known issues like this here
                    if ind==329 && trigs(eegInd)==160
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd).type = []; % delete
                    % update these
                     EEG.event(eegInd).edftype = []; %
                    trigs(eegInd) = [];
                     sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end

                     elseif regexp(fileNames{iB}, 'P15_2_2') % can fix known issues like this here
                    if ind==519 && trigs(eegInd)==58
                    % trig 13 missing, seems replaced by trg 45 (times
                    % match ET exactly, so overwriting it)
                    EEG.event(eegInd).type = 26; % delete
                    % update these
                     EEG.event(eegInd).edftype = 26; %
                    trigs(eegInd) = 26;
                    % sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end

                       elseif regexp(fileNames{iB}, 'P16_1_5') % can fix known issues like this here
                    if ind==863 && trigs(eegInd)==31
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd).type = []; % delete
                    % update these
                     EEG.event(eegInd).edftype = []; %
                    trigs(eegInd) = [];
                     sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end

                     elseif regexp(fileNames{iB}, 'P14_1_3') % can fix known issues like this here
                    if ind==406 && trigs(eegInd)==33
                    % trig 13 missing, seems replaced by trg 45 (times
                    % match ET exactly, so overwriting it)
                    EEG.event(eegInd).type = []; % delete
                    % update these
                     EEG.event(eegInd).edftype = []; %
                    trigs(eegInd) = [];
                    sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end

                       elseif regexp(fileNames{iB}, 'P14_1_3') % can fix known issues like this here
                    if ind==732 && trigs(eegInd-1)==32
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd-1).type = []; % delete
                    % update these
                     EEG.event(eegInd-1).edftype = []; %
                    trigs(eegInd-1) = [];
                     sTimes(eegInd-1) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end

                      elseif regexp(fileNames{iB}, 'P16_2_4') % can fix known issues like this here
                    if ind==317 && trigs(eegInd)==32
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd).type = []; % delete
                    % update these
                     EEG.event(eegInd).edftype = []; %
                    trigs(eegInd) = [];
                     sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end


                     elseif regexp(fileNames{iB}, 'P16_3_3') % can fix known issues like this here
                    if ind==535 && trigs(eegInd)==32
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd).type = []; % delete
                    % update these
                     EEG.event(eegInd).edftype = []; %
                    trigs(eegInd) = [];
                     sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end


                     elseif regexp(fileNames{iB}, 'P16_3_4') % can fix known issues like this here
                    if ind==364 && trigs(eegInd)==31
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd).type = []; % delete
                    % update these
                     EEG.event(eegInd).edftype = []; %
                    trigs(eegInd) = [];
                     sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end

                    elseif regexp(fileNames{iB}, 'P15_3_5') % can fix known issues like this here
                    if ind==840 && trigs(eegInd)==40
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd).type = []; % delete
                    % update these
                     EEG.event(eegInd).edftype = []; %
                    trigs(eegInd) = [];
                     sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end


                    elseif regexp(fileNames{iB}, 'P17_2_2') % can fix known issues like this here
                    if ind==807 && trigs(eegInd)==58
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd).type = []; % delete
                    % update these
                     EEG.event(eegInd).edftype = []; %
                    trigs(eegInd) = [];
                     sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end

                      elseif regexp(fileNames{iB}, 'P17_2_3') % can fix known issues like this here
                    if ind==364 && trigs(eegInd)==47
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd).type = 15; % delete
                    % update these
                     EEG.event(eegInd).edftype = 15; %
                    trigs(eegInd) = 15;
                    % sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end

                       elseif regexp(fileNames{iB}, 'P17_2_3') % can fix known issues like this here
                    if ind==459 && trigs(eegInd)==35
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd).type = []; % delete
                    % update these
                     EEG.event(eegInd).edftype = []; %
                    trigs(eegInd) = [];
                     sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end

                     elseif regexp(fileNames{iB}, 'P17_2_3') % can fix known issues like this here
                    if ind==752 && trigs(eegInd-1)==47
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd-1).type = 15; % delete
                    % update these
                     EEG.event(eegInd-1).edftype = 15; %
                    trigs(eegInd-1) = 15;
                    % sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end
elseif regexp(fileNames{iB}, 'P17_2_6') % can fix known issues like this here
                    if ind==399 && trigs(eegInd)==32
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd).type = []; % delete
                    % update these
                     EEG.event(eegInd).edftype = []; %
                    trigs(eegInd) = [];
                     sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end


                    elseif regexp(fileNames{iB}, 'P17_2_6') % can fix known issues like this here
                    if ind==591 && trigs(eegInd-1)==190
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd-1).type = []; % delete
                    % update these
                     EEG.event(eegInd-1).edftype = []; %
                    trigs(eegInd-1) = [];
                     sTimes(eegInd-1) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end

                    elseif regexp(fileNames{iB}, 'P18_1_1') % can fix known issues like this here
                    if ind==832 && trigs(eegInd)==35
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd).type = []; % delete
                    % update these
                     EEG.event(eegInd).edftype = []; %
                    trigs(eegInd) = [];
                     sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end

                      elseif regexp(fileNames{iB}, 'P14_2_3') % can fix known issues like this here
                    if ind==380 && trigs(eegInd)==48
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd).type = []; % delete
                    % update these
                     EEG.event(eegInd).edftype = []; %
                    trigs(eegInd) = [];
                     sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end

                    elseif regexp(fileNames{iB}, 'P14_2_3') % can fix known issues like this here
                    if ind==406 && trigs(eegInd-1)==59
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd-1).type = 27; % delete
                    % update these
                     EEG.event(eegInd-1).edftype = 27; %
                    trigs(eegInd-1) = 27;
                    % sTimes(eegInd-1) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end

                      elseif regexp(fileNames{iB}, 'P14_3_3') % can fix known issues like this here
                    if ind==628 && trigs(eegInd)==35
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd).type = 3; % delete
                    % update these
                     EEG.event(eegInd).edftype = 3; %
                    trigs(eegInd) = 3;
                    % sTimes(eegInd-1) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end

                    elseif regexp(fileNames{iB}, 'P14_2_4') % can fix known issues like this here
                    if ind==173 && trigs(eegInd)==41
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd).type = []; % delete
                    % update these
                     EEG.event(eegInd).edftype = []; %
                    trigs(eegInd) = [];
                     sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end

                    elseif regexp(fileNames{iB}, 'P14_2_6') % can fix known issues like this here
                    if ind==188 && trigs(eegInd)==45
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd).type = []; % delete
                    % update these
                     EEG.event(eegInd).edftype = []; %
                    trigs(eegInd) = [];
                     sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end


                      elseif regexp(fileNames{iB}, 'P17_3_5') % can fix known issues like this here
                    if ind==556 && trigs(eegInd)==175
                    % there are 121 conf-resps recorded, but only 120 in mat+edf. 
                    % it was trial 56, there were 2 confidence-responses
                    % recorded, a '8' 1ms after the 'n', and then a '9' ~350ms
                    % later. so I am going to delete that first one, to align
                    % it with mat+edf.
                    EEG.event(eegInd).type = []; % delete
                    % update these
                     EEG.event(eegInd).edftype = []; %
                    trigs(eegInd) = [];
                     sTimes(eegInd) = [];

                    % check issue is fixed
                    eegRelTrig = ismember(trigs, triggers.relevant);
                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
                    if any(tr ~= tr2) || any(freq ~= freq2)
                        disp([tr, freq, freq2, freq~=freq2])
                        keyboard('still a trigger mismatch');
                    end
                    end
 % 
 % elseif regexp(fileNames{iB}, 'P10_3_4') % can fix known issues like this here
 %                    if ind==154 && trigs(eegInd-1)==45
 %                    % trig 13 missing, seems replaced by trg 45 (times
 %                    % match ET exactly, so overwriting it)
 %                    EEG.event(eegInd-1).type = 26; % delete
 %                    % update these
 %                     EEG.event(eegInd-1).edftype = 26; %
 %                    trigs(eegInd-1) = 26;
 %                    % sTimes(eegInd) = [];
 % 
 %                    % check issue is fixed
 %                    eegRelTrig = ismember(trigs, triggers.relevant);
 %                    [freq, tr] = CountUnique(trigs(eegRelTrig)');
 %                    if any(tr ~= tr2) || any(freq ~= freq2)
 %                        disp([tr, freq, freq2, freq~=freq2])
 %                        keyboard('still a trigger mismatch');
 %                    end
 %                    end
                else
                    keyboard; % fix it
                end
    
            end
        end

        %% find epochs
                
        % find indices of triggers corresponding to the events we want to
        % time-lock to (evidence onset)
        trialStarts = find(trigs==triggers.trialStart);
        
        % there are some missing or misaligned triggers (i.e. mismatched to
        % trials), I have manually fixed these known issues here, and am
        % matching trigger RT to behavioural RT below to check trials are
        % correctly aligned
        if any(length(trialStarts) ~= [length(mat.resp.time) sum(~isnan(mat.resp.time))])
%             keyboard;
%                 error('wrong n trigs');
            trialStarts = find(trigs==triggers.trialStart);
        end
               

        % Now loop through the single trials and grow the trial-descriptor vectors and 'erp' matrix
        for iE = 1:length(trialStarts)
            
            nTrSoFar = sum(nTrPerBlock .* ((1:fileInfo.nBlocks)<iB)); % set prev blocks to zero
            iTr = nTrSoFar + iE; % trial index in entire experiment

            %store beh data now, so will still be there even if trial is
            %invalid later
            dataMat(iTr,:) = origMat(iE,:); 
            
            % first, check whether the epoch relative to the current event is even contained within the bounds of the EEG data (maybe the recording was stopped mid-trial?) - if not, skip the trial ('continue')
            if sTimes(trialStarts(iE)) < abs(eeg.epochSamples(1)) 
                continue; 
            end
            if (sTimes(trialStarts(iE)) + eeg.epochSamples(end)) > size(EEG.data,2)
%                 if regexp(fileNames{iB}, 'P01_3_1.bdf|P03_3_5.bdf'); continue; end
% %                 keyboard;
                continue;
%                 error('trigs errors');
            end
            
            thisTr = trialStarts(iE); % start of this trial

            % up until next trial/end of block
            if iE<length(trialStarts) % check in this epoch
                nextTr = trialStarts(iE+1);
            else
                nextTr = length(trigs); % if last, check until end
            end

            % find ev onset
            thisTarg = -1 + thisTr + find(ismember(trigs(thisTr:nextTr-1), triggers.evOnset),1);
            if isempty(thisTarg) % if no evOnset found, is early-resp
                if mat.resp.perf(2,iE)==mat.par.TOOEARLY
                    disp('response too-early, skipping trial');
                    continue;
                end
                keyboard;
%                 error('ev onset');
            elseif (sTimes(thisTarg) + eeg.epochSamples(end)) > size(EEG.data,2)
                % end of trial not within data
%                 keyboard;
                continue; % skip to next trial
            end
            

            % can be response or stim off
            nextRespInd = -1 + thisTarg + find(ismember(trigs(thisTarg:nextTr-1), triggers.resp),1);

            % get confResp here, to check response is before it
            stimOffInd = find(ismember(trigs, triggers.stimOffset) & sTimes>sTimes(thisTarg) & sTimes<=sTimes(nextTr),1);

            % resp is in this trial, and a non-NaN RT was recorded
            if ~isempty(nextRespInd) && nextRespInd<length(trigs) && ismember(trigs(nextRespInd), triggers.resp) && ~isnan(mat.resp.time(iE)) && (isempty(stimOffInd) || nextRespInd < stimOffInd)
                RS(iTr) = sTimes(nextRespInd)-sTimes(thisTarg);
                initResp(iTr) = find(trigs(nextRespInd)==triggers.resp); % 1-6 button
                certainty(iTr) = abs(round(initResp(iTr) -3.5)); % 1-6 -> [3 2 1 1 2 3
                respLR(iTr) = (initResp(iTr) > 3) + 1; % l/r
            else
                continue; % finish here
                % should store triggers anyway?
            end    
            
            % check RT roughly matches respT
            if abs( (RS(iTr)/eeg.fs*1000) - mat.resp.time(iE)*1000 ) > 15 % ms
                keyboard;
%                 error('RT mismatch');
            end
                        
                        
            % get trial info
            blockNum(iTr) = fileInfo.blockNums(thisPP(iB));
    
            % get stimOnset trig, to find corrLR
            stimOnTrigInd = -1 + thisTr + find(ismember(trigs(thisTr:thisTarg), triggers.evOnset),1);
            corrLR(iTr) = 2 - any(trigs(stimOnTrigInd) == triggers.isLeft); %1=left, 2=right

            % find cue trigger
            cueTrigInd = -1 + thisTr + find(ismember(trigs(thisTr:thisTarg), triggers.cueOff),1);
            cues(iTr) = find(trigs(cueTrigInd) == triggers.cueOff); % 1=neutral, 2=valid, 3=invalid


            %% Now extract the epoch
            ep = EEG.data(:,sTimes(thisTarg) + eeg.epochSamples);
            %Baseline correction
            blAmp(:,1,iTr) = nanmean(ep(:,isBetween(eeg.epochTimes, eeg.blWin)),2); % store this
            ep = ep - blAmp(:,1,iTr);
            erp(:,:,iTr) = ep; % now add the current epoch onto the growing 'erp' matrix by concatenation along the 3rd dimension

            %%% uncomment this if doing conf-resp locking
%             %% only keep the bits around the confidence response times
%             % [-1500 500];
%             
%             % check there was a confresp
%             if ~isempty(confCueInd) && ~isempty(confRespInd)
% 
%                 if (sTimes(confRespInd) + eeg.confRespSamples(end)) <= size(EEG.data,2)
%                     ep = EEG.data(:, sTimes(confRespInd) + eeg.confRespSamples);
%                 else
%                     inds = eeg.confRespSamples <= (size(EEG.data,2) - sTimes(confRespInd));
%                     ep = EEG.data(:, sTimes(confRespInd) + eeg.confRespSamples(inds));
%                 end
%                 ep = ep - blAmp(:,1,iTr); % baseline using the ev-locked baseline
% 
%             else
%                 ep = NaN(eeg.nChansTot, length(eeg.confRespSamples));
%             end
% 
%             erp(:,1:size(ep,2),iTr) = ep; % now add the current epoch onto the growing 'erp' matrix by concatenation along the 3rd dimension
            
            %% also store al the triggers in an epoch
            
            %%% uncomment this if doing confResp locking
%             if ~isempty(confCueInd) && ~isempty(confRespInd) 
%                 % get all sTimes within the epoch
%                 trigInds = isBetween(sTimes, sTimes(confRespInd) + minMax(eeg.confRespSamples,2));

            trigInds = isBetween(sTimes, sTimes(thisTarg) + minMax(eeg.epochSamples,2));
            epochTrigs{iTr} = trigs(trigInds);
            epochSTimes{iTr} = sTimes(trigInds) - sTimes(thisTarg);
%             end
            
            %% store saccade info
            saccAmpls{iTr} = [EEG.event(trigInds).sac_amplitude]; 
            saccAmpls{iTr}(~ismember(epochTrigs{iTr}, [triggers.eyes{1,2} triggers.eyes2{1:2,2}])) = NaN; % set non-saccade events to NaN


        end
    end
    
    disp(iTr);
    
    %% It's worth taking a look at this subject's VEOG, just for the last block, to get a sense whether blinks will be well
%     % detected by a certain threshold:

    isLastBlock = blockNum==max(blockNum);
    clf;
    plot(reshape(erp(129,:,isLastBlock), 1,[]));
    yline(250,'-k');
    title(['VEOG ' fileInfo.ppID{iPP}])
    
    %% Now that we have all trials from all blocks, save everything for this subject:

    % save EEG data separately, as can replace just that with interp/csd
    epochs.erp = erp(:,:,1:iTr);
    epochs.baseline = blAmp(:,:,1:iTr);
    epochs.blockNum = blockNum(1:iTr); % duplicate this

    
    % save all behavioural data separately then, to avoid repetition
    beh.blockNum = blockNum(1:iTr);
    beh.isLeft = corrLR(1:iTr);
    beh.cues = cues(1:iTr);
    beh.respLR = respLR(1:iTr);
    beh.RS = RS(1:iTr); % maybe convert into ms from samples?
    beh.corrLR = corrLR(1:iTr);
    beh.initAcc = double(respLR(1:iTr) == corrLR(1:iTr));
    beh.initResp = initResp(1:iTr);
    beh.certainty = certainty(1:iTr);
    beh.dataMat = dataMat(1:iTr,:);
    beh.dataMatNames = dataMatNames;
    beh.saccAmpls = saccAmpls(1:iTr);
    beh.trigs = epochTrigs(1:iTr);
    beh.sTimes = epochSTimes(1:iTr);
    
%     %% append only dataMat
%     AppendIntoFiles(fileInfo.ppID{iPP}, epochs)
     
    %% save EEG data
    fprintf('\n Saving raw file...\n');
    save(fullfile(fileInfo.rawFolder, [fileInfo.ppID{iPP} '_raw']),...
        '-v7.3', '-struct','epochs');

    %% save behavioural + triggers (just not EEG)

    fprintf('\n Saving beh file...\n');
    save(fullfile(fileInfo.behFolder, [fileInfo.ppID{iPP} '_beh']),...
        '-v7.3', '-struct','beh');

fprintf('\n end\n');

end