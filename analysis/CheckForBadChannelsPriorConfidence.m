% CheckForBadChannelsPriorConfidence
% This script opens matfiles with single-trial ERPs that have been
% extracted - usually without any re-referencing or interpolation ('raw')
% but possibly with low-pass filter and detrending done - and examines the
% data for channels that are abnormally high in variance
clear; close all;

outFolder = 'D:\cue_task\analysis\Data\Saves';
load(fullfile(outFolder, 'ExtractEpochsPriorConfidence.mat'));

interpPerBlock = 1; % do the checks per block
if interpPerBlock
    loadName = 'ChansToInterpBlock.mat';
else
    loadName = 'ChansToInterp2.mat';
end

loadOld = 0;
if loadOld && exist(loadName, 'file')
    % load this - will have ones already done
    c = load(loadName, 'chansToInterp');
    chansToInterp = c.chansToInterp; 
    
    % don't do the ones already done?
    ppsToDo = find(~ismember(fileInfo.ppID, fieldnames(c.chansToInterp)))';
%     ppsToDo = 1:fileInfo.nPP;
else
    chansToInterp = struct();
    ppsToDo = 1%:fileInfo.nPP;
end



%% 

% We're going to plot the standard deviation of each channel and spot ones that stick out. It's worthwhile to particularly
% scrutinize electrodes that might be particularly important in the final analysis. Here, we want to use electrodes around
% electrode CPz and around left and right motor areas (see 'cap_128_layout_medium.jpg'):
importantChans = [2 3 4 19 84 85 86 52 53 54 110 114 115 104 103 94 93 81 80 72 71 59]; % so we can pay particular attention to those
reRefChan = 23; % pick a channel to re-reference the data to and get a second picture of variance across channels - 
% important just because sometimes the reference used to load the EEG might itself have been bad during the recording...
% This is Oz

% Cz=A1, Pz=A19, Oz=A23, FPz=C17, Fz=C21, C3=D19, T7=D23, C4=B22, T8=B26
% [1, 19, 23, 54, 58, 81, 85, 115, 119]
%%

% Manually go through each subject in turn and use the code snippets that are commented out near the end to investigate bad channels
for iPP = 1 % change the number manually, before the %, and hit F5
    close all;
    disp(fileInfo.ppID{iPP});
    data = load(fullfile(fileInfo.rawFolder, [fileInfo.ppID{iPP} '_raw']),'erp','blockNum'); % contains epochs structure's fields
    
    data.erp = data.erp(1:eeg.nChans,:,:); % only EEG channels

    %%
    blockNums = unique(data.blockNum(~isnan(data.blockNum)));
    
    % concatenate all epochs and get Standard deviation (SD) per channel
%     [SD, SD2] = deal(zeros(eeg.nChansTot, length(unique(blockNum))));
%     for iBl = blockNums
%         iTrs = find(blockNum==iBl); % find the indices of the trials from the current block, which are not outliers
%         conc = reshape(erp(:,:,iTrs),[eeg.nChansTot,length(iTrs)*eeg.nSamples]); % concatenate all trials
%         conc2 = conc - repmat(conc(reRefChan,:),[eeg.nChansTot,1]); % Also look at SD when referenced to somewhere at the back - electrode Oz for example
%         for iCh = 1:eeg.nChansTot % include externals here as well, because it might be a good idea to check whether those are noisy too
%             SD(iCh,iBl) = std(conc(iCh,:)); % measure S.D. of each channel 
%             SD2(iCh,iBl) = std(conc2(iCh,:));
%         end
%     end

    %% remove very bad trials
    
    % % %     % can remove entire bad trials and recalc SD?
    badTrials = find(nanmean(nanstd(data.erp(1:eeg.nChans,:,:),[],2) > 50, 1)>.7);
    if any(badTrials)
        disp(badTrials);
        data.erp(:,:,badTrials) = NaN;
    end

    % split by block - not groupmeans in case just one block
    blocks = unique(data.blockNum);
    blocks(isnan(blocks)) = [];
    blocks(blocks==0) = []; % ignore - not sure why these are
%     blocks = 1:6;
    nBlocks = max(blocks); % ignore any where block==0
    conc = cell(1,nBlocks);
    [SD, SDr] = deal(NaN(eeg.nChans, max(nBlocks)));
    for i = 1:max(nBlocks)
        x = reshape(data.erp(:,:,data.blockNum==blocks(i)),eeg.nChans,[]);
        SD(:,i) = nanstd(x,[],2);
        SDr(:,i) = nanstd(x - x(reRefChan,:),[],2);

%         conc{i} = data.erp(:,:,data.blockNum==blocks(i));
    end
%     conc = nancat(conc); %[chans time tr blocks]
%     nTr = size(data.erp,3);
%     conc = reshape(conc, eeg.nChans, eeg.nSamples*nTr, nBlocks);
%     conc = reshape(data.erp, eeg.nChans, eeg.nSamples*fileInfo.nTr, fileInfo.nBlocks); % [chan time*tr block]
    
%     SD = sq(nanstd(conc, [], 2));
%     SDr = sq(nanstd(conc - conc(reRefChan,:,:),[],2));
%     SDr = sq(nanstd(conc - nanmean(conc,1),[],2));
    
    
    %% also get artefact counts
    
    isArt = sq(max(abs(data.erp(1:eeg.nChans, isBetween(eeg.epochTimes, [-1400 1800]),:)),[],2) - 100 > 0);  % is uV over 100 in window?

    chanArts = sum(isArt,2); % arts per channel
    chanArtsPerBlock = NaN(eeg.nChans, nBlocks);
    for i = 1:nBlocks
        chanArtsPerBlock(:,i) = sum(isArt(:,data.blockNum==blocks(i)),2);
    end
%     chanArtsPerBlock = sq(sum(reshape(isArt,eeg.nChans,80,16),2));
    % will interpolate if > 40? 30?
    
    %% are there any channels that stick out in terms of variance?
    
    
    figure();
    plot(0:100, prctile(col(SD), 0:100), '-x');
    xlabel('%'); ylabel('SD');
    ylim([0 200]);
    % find the inflection point, can use that to set sdLims
    sdLims = [1 50]  
    
%     keyboard; % adjust sdLims
    %%
    figure; 
    subplot(2,1,1); hold on; plot(SD); 
    ylabel('SD'); xlabel('channel'); 
    title(['pp ' num2str(iPP) ' ' fileInfo.ppID{iPP}])
    plot([1; 1]*importantChans, [0;1] * max(SD(importantChans,:),[],2)', '-k')
    legend(num2str(blockNums'),'Location','Best');
    yline(sdLims(2),':k');
    ylim([0 200]);
    
    subplot(2,1,2); 
    hold on; 
    plot(SDr); 
    ylabel('SD'); 
    xlabel('channel');
    legend(num2str(blockNums'),'Location','Best');
    yline(sdLims(2),':k');
    ylim([0 200]);

    %%
    figure();
    imagesc(SD, [0 100] .* [1 1]); colorbar
    ylabel('channels'); xlabel('block');
    % horizontal stripes are bad channels, vertical are bad blocks
    
    SDct = squeeze(nanstd(data.erp,[],2));
    figure; imagesc(SDct,[0 100].*[1 1]); ylabel('channel'); xlabel('trial'); colorbar
    xticks(0:fileInfo.nTr:fileInfo.maxTr)

%     hold on; plot(cumsum(CountUnique(data.blockNum)),eeg.nChans, 'kx')
%     horizontal stripes are bad channels, vertical are bad trials/times
%     
%     figure(); 
%     imagesc(SDct(importantChans,:), sdLims.*[1 1.2]); 
%     hold on; plot(80:80:1280,max(ylim), 'kx')
%     yticks(1:length(importantChans)); yticklabels(importantChans); 
%     xticks(1:fileInfo.nTr:fileInfo.maxTr); xticklabels(1:fileInfo.nBlocks);
%     ylabel('important channels'); xlabel('block'); 
%     now taking a look through the marked 'important' channels, I see 52 sort of sticks out in block (check [sortedSD,I] = sort(SD(52,:)))...5! So finally: 
%     
%%     plot artefacts too
    figure();
    imagesc(chanArtsPerBlock, [0 240]); colorbar
    title('# artefacts');
    ylabel('channel'); xlabel('block');
    
%     [bl,ch] = find(chanArtsPerBlock' > 40); sortrows([bl,ch])


    % other useful plotting code:
    % plot a selected bunch of electrodes, concatenated across trials:
%     elec2plot = [2 14 17 23 36 55 56]; figure; hold on; for q=1:length(elec2plot), plot((q-1)*100+conc(elec2plot(q),:)); end; legend(num2str(elec2plot'))

    %% examine the plots
    
%     ch2interp = cell(1,length(blockNums)); % makes an empty cell array for filling in the bad channels for this subject
   
    % look at the figures - the first figure (subplot) can be used to see
    % which blocks stick out - the second one (imagesc) also shows this
    
    % the next one (imagesc) shows SD on each trial, look for horizontal or
    % vertical stripes to decide if it is a bad channel or trial
    % you can click on a bright patch and use (get, gca, 'CurrentPoint') to
    % find the [x y z] coords - there will be some slight rounding errors
    % depending on how large the image is onscreen
    
    % the next one shows the same but just for the 'importantChans'
    
    % use 
    % [bl,ch] = find(SD(1:128,:)'>sdLims(2)); sortrows([bl,ch])
    % to find the blocks that are bad for a certain channel
    
    % or
    % [b, i] = sort(SD(13,:))
    
    % add them to the 'bad' variable below = {block [badChansInThatBlock]}
    % This will fill in chansToInterp, which is saved for the
    % InterpBadChannels.m file
    
    


    [bl, ch] = find(~isBetween(SD', sdLims) | chanArtsPerBlock' > 111);
    sortrows([bl, ch]);

%      nBadChans = length(unique(ch))% > 10%?
%     [a,b] = CountUnique(bl); [b,a]

    ch2interpAll = [unique(ch)']  % all ones being done
    disp(eeg.chanNames(ch2interpAll)');


    if interpPerBlock 
        ch2interp = cell(1, nBlocks);
        for i = 1:nBlocks
            ch2interp{i} = ch(bl==i);

            fprintf('\n\nblock %d = %d chans:\n',i, length(ch2interp{i}));
            disp(eeg.chanNames(ch2interp{i})');
        end
        nPerBlock = cellfun(@length, ch2interp);
        if any(nPerBlock > 12)
            warning('>12 channels found per block: ')
            disp(nPerBlock)
        end
    else
        ch2interp = ch2interpAll;
    end

    
%     % use this to insert channels into all blocks
%     insertChans = 77;%ch2interpAll';%    [103];
%     for i=1:nBlocks
%         ch2interp{i} = unique([ch2interp{i}; insertChans]);
%     end
    disp(sq(nancat(ch2interp)));

%     
    %% show topoplot of them

    if ~exist('topoplot','file'); eeglab nogui; end
    figure();
    topoplot(nanmean(SD,2),eeg.chanlocs,'electrodes','off',...
        'emarker2', {ch2interpAll, '.','k',10,1},'mapLimits',[0 100]);%[0 max(nanmean(SD,2))]);
    colorbar;

    % plot each block
    figure;
    for i=1:max(nBlocks)
        subplot(2,3,i);
        topoplot(SD(:,i),eeg.chanlocs,'electrodes','off',...
            'emarker2', {ch2interp{i}, '.','k',10,1},'mapLimits',[1 100]);%[0 max(nanmean(SD,2))]);
        colorbar;
        title(['block ' num2str(i)]);
    end
% 
%     figure
%     h = errorBarPlot(permute(data.erp(ch2interp,:,:),[3 2 1]),'area',1,'xaxisvalues',eeg.epochTimes); legend([h{:,1}], eeg.chanNames(ch2interp));
%     % uncomment below to show some good channels
%     h = errorBarPlot(permute(data.erp(find(~ismember(1:128,ch2interp),10),:,:),[3 2 1]),'area',1,'xaxisvalues',eeg.epochTimes); legend([h{:,1}], eeg.chanNames(find(~ismember(1:128,ch2interp),10)));

    keyboard; % edit the thresholds above, or add others into [bl ch] to be included below
%% for picking the baddest ch

% for i = 1: 6
% if length (ch2interp{i}) >12
% a = SD (ch2interp{i}, i);
% j = [ch2interp{i}, a];
% k = sortrows (j, -2);
% tosee = k (1:12, 1);
% tosee = sortrows (tosee, 1);
% ch2interp{i} = tosee;
% end
% end
% 
% keyboard
%% add chans to interpolate per block

%     for iBl = 1:fileInfo.nBlocks
%         nBadChansPerBlock = length(unique(ch(bl == iBl)));
%         if nBadChansPerBlock > eeg.nChans/10
%             keyboard;
%         end
%         ch2interp{iBl} = ch(bl == iBl);
%     end
    
    %% store
   
    chansToInterp.(fileInfo.ppID{iPP}) = ch2interp; % store each pp in a field, cell per block if doing that   
    chansToInterp = orderfields(chansToInterp);
    
    save(loadName, 'chansToInterp', 'importantChans');

    %%%% also store in the loaded file
    save(fullfile(fileInfo.rawFolder, [fileInfo.ppID{iPP} '_raw']),...
        '-append','ch2interp','sdLims');
end


