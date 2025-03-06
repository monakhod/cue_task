% DoFFTDelayedConf
% run per person, just get beta for now, as can use smaller window then
% save per person. Instead of averaging, just save all freqs <50hz
% for useCSD = 1 %1:2 will do both

useCSD = 1;
clc; close all; clearvars -except useCSD;
outFolder = 'D:\cue_task\analysis\Data\Saves';
outFolder2 = 'D:\cue_task\analysis\Data\cpp';

load(fullfile(outFolder, 'ExtractEpochsPriorConfidence.mat'));
if useCSD
    fileInfo.fftFolder = 'D:\cue_task\analysis\Data/FFTCSD/';
    fileInfo.interpFolder = 'D:\cue_task\analysis\Data\General';
    saveName1 = 'DoFFTDelayedConfCSD.mat';
    sortName = '_cpps.mat';
    intpName = '_whole.mat';
end
load(fullfile(outFolder, 'BehDataLoadPriorConfidence.mat'),'behDataSess'); % for RS

% remove artefacts
% load(fullfile(outFolder, 'FlagArtefactsDelayedConf_Long.mat'),'isFlagged');


%% SET
setID = {'P22'};
for p = 1: length (setID)
fileInfo.ppID = {setID{p}};
fileInfo.maxTr = 2592;

fs=eeg.fs;
% STFT parameters
nFreqs = round( (fs/25) * 10); % length of window = num freqs,
freqs = ((0:nFreqs-1)*fs) / nFreqs; % frequency scale, given window length (remember resolution = 1/window-duration)
windowSize = nFreqs; % samples; 500ms; 4 cycles of 8Hz

freqBands = [0 40]; % will save all within this


% stim windows, must be windowSize from edges of epoch
stimWindows = -1400:50:1800; % in msec, centered on what times do you want to measure spectral amplitude? i.e., where to center each consecutive window in time

% resplocked
respWindows = -1400:50:50; % so take window 500ms either side
% these are for the resp-locking - i.e. the erp
respWindowS = 1639 + [-1638:256]; % eeg.epochTimes starts at -1600, so offset by that
respTimes = ((-1638:256) ./ 1024) .* 1000; % times
nRespTimes = length(respWindowS);


save(fullfile(outFolder,saveName1), 'fs','nFreqs','freqs','freqBands','windowSize',...
    'stimWindows','respTimes','respWindows')

%%
for iPP = 1:fileInfo.nPP
if 1%~exist( fullfile(fileInfo.fftFolder, [fileInfo.ppID{iPP} '_resp_conf.mat']), 'file')

    disp(iPP);
    tic;
    data = load(fullfile(fileInfo.interpFolder, [fileInfo.ppID{iPP} intpName]),'erp_gen'); % interpolated, ev-locked epochs
    % erp = single
load(fullfile(outFolder2, [fileInfo.ppID{iPP} sortName]),'RT'); % interpolated, ev-locked epochs
%     data.erp(:,:,isFlagged(:,iPP)) = NaN; % remove flagged

    %% compute short time fourier transform (STFT) for each single trial:
    STFT  = myFFT(data.erp_gen(1:eeg.nChansTot,:,:), stimWindows, eeg.epochTimes, freqs, freqBands);

    % save as single?
    STFT = single(STFT);

    save(fullfile(fileInfo.fftFolder, [fileInfo.ppID{iPP} '_stim']),'-v7.3',...
        'STFT','freqBands');
        

    %% also do resp locked
    % quicker to recalc this than load it up

    respErp = NaN(eeg.nChansTot, length(respTimes), fileInfo.maxTr);
    % get from -1250:1250 around response
    for i = 1:fileInfo.maxTr
        if ~isnan(RT(i))% && ~isFlagged(i, iPP)
            windowInds = RT(i) + respWindowS;
            if all(isBetween(minMax(windowInds,2), [1 eeg.nSamples]))
                respErp(:,:,i) = data.erp_gen(1:eeg.nChansTot, windowInds, i); % whole epoch is used
            else
                % if the RT is < 250ms then will have to nanpad?
                windowInds = windowInds(isBetween(windowInds,[1 eeg.nSamples])); % keep valid
                respErp(:,:,i) = [NaN(eeg.nChansTot,length(respTimes)-length(windowInds))  data.erp_gen(1:eeg.nChansTot,windowInds,i)];
            end
        end
    end
    % run SFTF on this
    STFT = myFFT(respErp, respWindows, respTimes, freqs, freqBands);
    
    % save as single?
    STFT = single(STFT);
    
    save(fullfile(fileInfo.fftFolder, [fileInfo.ppID{iPP} '_resp']),'-v7.3',...
        'STFT','freqBands');
    
    t = toc

   

end    
end

% end
end
function STFT  = myFFT(erp, windowCentres, epochTimes, freqs, freqBands)
% STFT = myFFT(erp, windowCentres, epochTimes, freqs, freqBands)
% do FFT on each time window, get abs power, take mean within frequency
% bands
% Inputs:
%   erp = [nChans, nTimes, nTrials]
%   windowCentres = time points of centres of windows
%   epochTimes = time points of erp [1 nTimes]
%   freqs = vector of frequencies desired, will be used to pick window
%       widths (1 sample per frequency)
%   freqBands = [min max] of all freqs to save
% 
% Outputs:
%   STFT = [nChans, nWindows, nTrials, nFreqsInBand] mean abs power in each
%       freq within band


nFreqs = length(freqs);
bandInds = isBetween(freqs, freqBands);
nB = sum(bandInds);

[nChans,~,nTr] = size(erp);
STFT = NaN(nChans, length(windowCentres), nTr, nB); % preallocate

for iT = 1:length(windowCentres) % for each STFT timepoint (window centre) we'll compute the FFT for all trials at once
    
    
    [~,samp] = min(abs(epochTimes-windowCentres(iT))); % find the sample point in the ERP epoch corresponding to the centre of the current FFT window
    spec = abs( fft( erp(:,samp-round(nFreqs/2)+(1:nFreqs),:) ,[],2) ) ./ (nFreqs/2); % compute the magnitude of FFT (and scale it so that it's amplitude in same units as EEG
    %[ch freqs tr]
    STFT(:,iT,:,:) = permute(spec(:,bandInds,:),[1 3 2]);
    
end

end
