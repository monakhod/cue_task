% GetMuBetaDelayedConf
% load FFT - STFT 4th dim is freq band, [beta alpha/mu theta] - can combine
% beta & alpha/mu if want 8-30, otherwise beta is just 14-30hz
% remove flagged, grand topo pre-resp, pick elecs, outlier removal
% go back and apply the same to the stimlocked too
clc; clear all; 

outFolder = 'D:\cue_task\analysis\Data\Saves';
outFolder2 = 'D:\cue_task\analysis\Data\General';
% outFolder = '../../../../../OneDrive - TCDUD.onmicrosoft.com/ProjectsOD/DelayedConf/Analysis/Saves/';



load(fullfile(outFolder, 'ExtractEpochsPriorConfidence.mat'))

pp = {'P01', 'P02','P03','P04','P05','P06','P07','P08','P10','P13',...
   'P14','P15','P16','P17','P18', 'P19', 'P20', 'P21', 'P22'};

for i = 1 :length (pp)

fileInfo.ppID = {pp{i}};
intpName = '_whole.mat';
sortName = '_sorted_t.mat';

load(fullfile(outFolder2, [fileInfo.ppID{1}, intpName]),'isGood_gen');

load(fullfile(outFolder2,  [fileInfo.ppID{1}, sortName]),'lr_resp'); % for l/r

loadUp = 0; % load up PreRespMeanVoltage file?
useCSD = 1; % CSD or voltage

if useCSD
%     fileInfo.fftFolder = '.\Data\FFTCSD';
    fileInfo.fftFolder = 'D:\cue_task\analysis\Data\FFTCSD';
    load(fullfile(outFolder, 'DoFFTDelayedConfCSD.mat')); % get FFT info
    loadName = 'PreRespMeanBetaCSD.mat';
    saveName = 'GetMuBetaDelayedConfCSD200.mat';
    mapLims = [-1.3 1.3];
end

freqs1 = freqs(isBetween(freqs, freqBands)); % trim to saved freqs

betaInds = find(isBetween(freqs1, [8 30]));
% betaInds(round(freqs1(betaInds))==20) = []; % remove 20hz (i.e. ssmep)

preRespWindow = [-150 -50];
preRespInds = isBetween(respWindows, preRespWindow);
fileInfo.maxTr = 2592;

%% load up all - get the grand average for good trials

if loadUp && exist(fullfile(outFolder, loadName), 'file')
    r = load(fullfile(outFolder, loadName),'betaLat');
    betaLat = r.betaLat;

else

    betaLat = NaN(eeg.nChans, fileInfo.nPP, 2);
    betaTopo = NaN(eeg.nChans, fileInfo.nPP, 2);
    for iPP = 1:fileInfo.nPP

        disp([fileInfo.ppID{iPP} '...'])

        % get resplocked data
        data = load(fullfile(fileInfo.fftFolder, [fileInfo.ppID{iPP} '_resp']),'STFT');

        % remove flagged trials
        data.STFT(:,:,isGood_gen==0,:) = NaN;
        
        % lateralisation = left resp, right resp
        betas = nanmean(nanmean(data.STFT(1:eeg.nChans,preRespInds,:,betaInds),4),2); % average over time + freq
        betaLat(:,iPP,:) = [nanmean(betas(:,1,lr_resp==1),3), nanmean(betas(:,1,lr_resp==2),3)];
        
        % also store this for lateralisation data for topoplot at RT
        betaTopo(:,iPP,:) = [nanmean(nanmean(data.STFT(1:eeg.nChans, respWindows==0, lr_resp==1, betaInds),4),3),...
            nanmean(nanmean(data.STFT(1:eeg.nChans, respWindows==0, lr_resp==2, betaInds),4),3)];

    end

    save(fullfile(outFolder, loadName), 'betaLat', 'preRespWindow','betaTopo');
end

%% combine across sessions first

% betaLat2 = permute(groupMeans(betaLat, 2, fileInfo.ppID1n,'dim'),[4,3,2,1]); %[ch pp resp sess]
% 
% % need to weight them by trials
% nTrs = permute(groupMeans(sum(~isFlagged), 2, fileInfo.ppID1n, 'dim'), [3,2,4,1]); %[1 pp1 1 sess]
% nTrs = repmat(nTrs, eeg.nChans, 1, 2, 1); %[ch pp 2 sess]
% 
% betaLat2 = nansum(betaLat2 .* nTrs, 4) ./ nansum(nTrs,4); % weighted average by trials
% 
% fileInfo.nPP1 = length(unique(fileInfo.ppID1n));


%% topoplot that

% Corbett et al., 2021, bioRxiv uses C3 + 4 (D19 + B22)
% plot average topography from -150:-50ms before response
% cppChanNames = {'D19','D18','D12'; 'B22','B21','B31'}'; % beta chans
% cppChanNames = {'D2','D13','C24','D12','D3'; 'C11','C3','C12','B32','B31'}'; % 

cppChanNames = {'D18','D19','D28', 'B22','B21', 'B18'}'; % 
% cppChanNames = {'D18','D14','D13', 'B21','B20','B32'}'; % 

[~, cppChanInds] = ismember(cppChanNames, eeg.chanNames);

% labels.response = {'left','right'};
% 
% if ~exist('topoplot','file')
%     eeglab nogui;
% end
% figure();
% t=tiledlayout('flow');ax=[];
% for i = 1:2
%     ax(i) = nexttile(t);
%     topoplot(nanmean(betaLat(:,:,i),2),...
%         eeg.chanlocs, 'electrodes','on','colormap',crameri('vik'),...
%         'emarker',{'.','k',10,1},'emarker2',{col(cppChanInds), '.','y',10,1},...
%         'mapLimits', [0 mapLims(2)*10]);
%     colorbar
%     title(labels.response{i});
% end
% 
% ax(3) = nexttile(t);
% topoplot(nanmean(diff(betaLat,[],3),2),...
%     eeg.chanlocs, 'electrodes','on','colormap',crameri('vik'),...
%     'emarker',{'.','k',10,1},'emarker2',{col(cppChanInds), '.','y',10,1});
% colorbar
% title('difference for ', fileInfo.ppID{iPP});
% % % pick two clusters (left + right hemisphere)
 chanClusters = reshape(cppChanNames,3,2);
 [~, chansInds] = ismember(col(chanClusters), eeg.chanNames);
% 
% difference = nanmean(diff(betaLat,[],3),2);
% 
% save(fullfile(outFolder2, [fileInfo.ppID{iPP}, '_difference_map_8to30']), 'cppChanNames', 'difference');
%% topoplot per person
% % 
% figure();
% t = tiledlayout('flow');
% for iPP = 1:fileInfo.nPP1
%     if 1%~ppsToExclude(iPP)
%         ax = nexttile(t);
%         topoplot(-diff(betaLat2(:,iPP,:),[],3),...
%             eeg.chanlocs, 'electrodes','off','colormap',crameri('vik'),...
%             'emarker',{'.','k',10,1}, 'emarker2', {chansInds, '.','k',10,1});
%         colorbar;
%         title(iPP);
%     end
% end
% keyboard; %%%%%% edit these picked channels

%% pick maximal in that per person

% betas = -diff(betaLat(chansInds,:,:),[],3); % get mean voltages per channel per pp   
% betas = reshape(betas, 3,2,[]); % left/ right
% 
% % pick maximal amplitude surrounding response execution    
% 
% flips = [1 -1]; % invert right side
% % do per side
% for i = 1:2
%     [m,j] = max(betas(:,i,:) * flips(i));
%     betasPicked(:,i)  = m *flips(i);
%     betaChans(:,i) = chanClusters(j,i); % get name
%     betaChanInds(:,i) = chansInds(j + (i-1)*3); % index
% end
%     

%% take mean in each cluster?

[~, clusterInds] = ismember((chanClusters), eeg.chanNames)

%% load that data again and store cpp chans, remove outliers

betas = struct();
betas.resp = NaN(fileInfo.nPP, length(respWindows), fileInfo.maxTr, 2); % 4th dim is left/right
betas.stim = NaN(fileInfo.nPP, length(stimWindows), fileInfo.maxTr, 2); % 4th dim is left/right
% betas.confResp = NaN(fileInfo.nPP, length(confRespWindows), fileInfo.maxTr, 2); % 4th dim is left/right

for iPP = 1:fileInfo.nPP
    disp([fileInfo.ppID{iPP} '...'])
    data = load(fullfile(fileInfo.fftFolder, [fileInfo.ppID{iPP} '_resp']),'STFT');
    nTrs = size(data.STFT,3);
    % store that channel as cpp
%     betas.resp(iPP, :, 1:nTrs,:) = permute(data.STFT(betaChanInds(iPP,:), :, :,betaInd),[4,2,3,1]);
    for i = 1:2 % take mean in each cluster
        betas.resp(iPP, :, 1:nTrs,i) = nanmean(nanmean(data.STFT(clusterInds(:,i), :, :, betaInds),4),1); % av over freqs, then chans
    end
    
    betas.resp(iPP,:,isGood_gen==0,:) = NaN; % remove flagged trials before mean calc
    
    
    %%%%%%% also do the stim-locked
%     disp([fileInfo.ppID{iPP} '...'])
    data = load(fullfile(fileInfo.fftFolder, [fileInfo.ppID{iPP} '_stim']),'STFT');
    nTrs = size(data.STFT,3);

    % store that channel as cpp
%     betas.stim(iPP, :, 1:nTrs,:) = permute(data.STFT(betaChanInds(iPP,:), :, :,betaInd),[4,2,3,1]);
    for i = 1:2
        betas.stim(iPP, :, 1:nTrs,i) = nanmean(nanmean(data.STFT(clusterInds(:,i), :, :, betaInds),4),1);
    end
    betas.stim(iPP,:,isGood_gen==0,:) = NaN; % remove flagged trials before mean calc
    

    %%%%%% confresplocked
%     data = load(fullfile(fileInfo.fftFolder, [fileInfo.ppID{iPP} '_resp_conf']),'STFT');
%     nTrs = size(data.STFT,3);
% 
%     % store that channel as cpp
% %     betas.confResp(iPP, :, 1:nTrs,:) = permute(data.STFT(betaChanInds(iPP,:), :, :,betaInd),[4,2,3,1]);
%     for i = 1:2
%         betas.confResp(iPP, :, 1:nTrs,i) = nanmean(nanmean(data.STFT(clusterInds(:,i), :, :, betaInds),4),1);
%     end
%     betas.confResp(iPP,:,isFlagged(1:nTrs,iPP),:) = NaN; % remove flagged trials before mean calc
    
end

%% make both into [ipsi contra ipsi-contra]

% flip those where respLR==2
for iPP = 1:fileInfo.nPP
    betas.resp(iPP,:,lr_resp==2,:) = flip(betas.resp(iPP,:,lr_resp==2,:),4);
    betas.stim(iPP,:,lr_resp==2,:) = flip(betas.stim(iPP,:,lr_resp==2,:),4);
    % betas.confResp(iPP,:,lr_resp==2,:) = flip(betas.confResp(iPP,:,lr_resp==2,:),4);
end
% 
% contra - ipsi
betas.resp(:,:,:,3) = diff(betas.resp,[],4);
betas.stim(:,:,:,3) = diff(betas.stim,[],4);
% betas.confResp(:,:,:,3) = diff(betas.confResp,[],4);

%% save

save(fullfile(outFolder,saveName), ...
    'betas','chanClusters','preRespWindow',...
    'preRespInds','respWindows','stimWindows');

save(fullfile(outFolder2,[fileInfo.ppID{iPP}, 'mu_beta']), ...
    'betas','chanClusters','preRespWindow',...
    'preRespInds','respWindows','stimWindows');
disp ('done');
end