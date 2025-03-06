% BasicPipelineStepsPriorConfidence
% dbstop if error;
% make folder names to process

p = 22;
pid = cellfun(@(x) sprintf('P%02d', x), num2cell(p)', 'Uni',0);

% pid = {'Debugging'};

%% check for mismatch
my_CheckForMismatch(pid);

%% convert edf2asc - all folders where it hasn't been run yet

edf2ascLoopPriorConfidence(pid);

%% Extraction

ExtractEpochsPriorConfidence(pid); % load up data, filter, epoch, save

%% find bad channels to mark for interpolation

CheckForBadChannelsPriorConfidence

%% interpolate - per block now

InterpBadChannelsBlock('D:\cue_task\analysis\Data\Saves\ExtractEpochsPriorConfidence.mat','D:\cue_task\analysis\Data\Interp');
% 
%% behavioural data

BehDataLoadPriorConfidence % doesn't remove flagged

%% Flag artefacts

FlagArtefactsPriorConfidence % loads beh data

%% CSD

CSDTransform('D:\cue_task\analysis\Data\Saves\ExtractEpochsPriorConfidence.mat'); %messing here, go to line  100

%% resplock
% see Kev's paper for his epochs - cue, stim, resp locked epochs
% can adapt flagging to run separately on each if desired

 %RespLockPriorConfidence
% 
% 
% 
% %% CPP
% % 
% GetCPPPriorConfidence
% GetCPPTopoPriorConfidence
% % 
% % % %% FFT
% % 
DoFFTDelayedConf
% % 
% % %% beta
% % 
GetMuBetaDelayedConf
% % % MuBetaAnalysisPriorConfidence
% % 
% % %% CNV?
% % 
% % %% SSMEP?
% % 
% % 
