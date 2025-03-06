clc; clear all; 

outFolder = 'D:\cue_task\analysis\Data\Saves';
outFolder2 = 'D:\cue_task\analysis\Data\General';
% outFolder = '../../../../../OneDrive - TCDUD.onmicrosoft.com/ProjectsOD/DelayedConf/Analysis/Saves/';



load(fullfile(outFolder, 'ExtractEpochsPriorConfidence.mat'))

 pp = {'P01', 'P02','P03','P04','P05','P06','P07','P08','P10','P13',...
   'P14','P15','P16','P17','P18', 'P19', 'P20', 'P21', 'P22'};
diff_for_all = []; twft =[]; twt = [];
for i = 1 :length (pp)

fileInfo.ppID = {pp{i}};
intpName = '_whole.mat';
sortName = '_sorted_t.mat';

load(fullfile(outFolder2, [fileInfo.ppID{1}, intpName]),'isGood_gen');

load(fullfile(outFolder2,  [fileInfo.ppID{1}, sortName]),'freq_ch', 'e_or_c', 'freq_correct'); % for correct fr

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
tweny_Hz = find ( (abs (freqs-20) <2.8));
twenFive_Hz = find ( (abs (freqs-25) <2.8));


preStimWindow = [-1000 -600];
preStimInds = isBetween(stimWindows, preStimWindow);
fileInfo.maxTr = 2592;
eeg.blWin = [-1300]; %-1400 -1250
 

%% load up all - get the grand average for good trials

    for iPP = 1:fileInfo.nPP

        disp([fileInfo.ppID{iPP} '...'])

        % get resplocked data
        data = load(fullfile(fileInfo.fftFolder, [fileInfo.ppID{iPP} '_stim']),'STFT');
        data2 = load(fullfile(fileInfo.fftFolder, [fileInfo.ppID{iPP} '_resp']),'STFT');
        % remove flagged trials
        data.STFT(:,:,isGood_gen==0,:) = NaN;
      
        
        % lateralisation = left resp, right resp
        tw = data.STFT(1:eeg.nChans,:,:, tweny_Hz(2))./ ...
            (data.STFT(1:eeg.nChans,:,:, tweny_Hz(1)) +data.STFT(1:eeg.nChans,:,:, tweny_Hz(3))); % average over time + freq
         blAmp = nanmean(tw(:,stimWindows==eeg.blWin,:),2); % store this %%%tw(:,isBetween(stimWindows, eeg.blWin),:)
tw = log (tw./ blAmp);
% figure;
%  h = errorBarPlot(sq(permute (tw(23,:,:), [3 1 2])), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);

twf = data.STFT(1:eeg.nChans,:,:, twenFive_Hz(2))./ ...
            (data.STFT(1:eeg.nChans,:,:, twenFive_Hz(1)) +data.STFT(1:eeg.nChans,:,:, twenFive_Hz(3))); % average over time + freq

blAmp = nanmean(twf(:, stimWindows==eeg.blWin,:),2); % store thistwf(:,isBetween(stimWindows, eeg.blWin), :)
twf = log (twf./ blAmp);

%% for rep locked
tw2 = data2.STFT(1:eeg.nChans,:,:, tweny_Hz(2))./ ...
            (data2.STFT(1:eeg.nChans,:,:, tweny_Hz(1)) +data2.STFT(1:eeg.nChans,:,:, tweny_Hz(3))); % average over time + freq
         blAmp = nanmean(tw2(:,stimWindows==eeg.blWin,:),2); % store this %%%tw(:,isBetween(stimWindows, eeg.blWin),:)
tw2 = log (tw2./ blAmp);
% figure;
%  h = errorBarPlot(sq(permute (tw(23,:,:), [3 1 2])), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);

twf2 = data2.STFT(1:eeg.nChans,:,:, twenFive_Hz(2))./ ...
            (data2.STFT(1:eeg.nChans,:,:, twenFive_Hz(1)) +data2.STFT(1:eeg.nChans,:,:, twenFive_Hz(3))); % average over time + freq

blAmp = nanmean(twf2(:, stimWindows==eeg.blWin,:),2); % store thistwf(:,isBetween(stimWindows, eeg.blWin), :)
twf2 = log (twf2./ blAmp);
% figure;
%  h = errorBarPlot(sq(permute (twf(23,:,:), [3 1 2])), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
    end


ssvep.twenty_stim = nanmean (tw(23:24,:,:), 1);
ssvep.twentyFive_stim = nanmean (twf(23:24,:,:), 1);
ssvep.twenty_resp = nanmean (tw2(23:24,:,:), 1);
ssvep.twentyFive_resp = nanmean (twf2(23:24,:,:), 1);

cppChanNames = {'A23','A24'}';
    [~, cppChanInds] = ismember(cppChanNames, eeg.chanNames);

    
save(fullfile(outFolder2,[fileInfo.ppID{iPP}, '_SSVEP']), ...
    'ssvep','preStimWindow',...
    'preStimInds','stimWindows');
end
% 
