% this code glueing together all bits together
% first flaggingArt and my_beh_check should be done


clear all
clc

ERPfolder = 'D:\cue_task\analysis\Data\CSD'; 
fold = 'D:\cue_task\analysis\Data\Combine';

PP = {'P22'}
f = what(fold);

for j = 1:length (PP)
wh = PP{j};
% files = f.mat(cellRegexpi(f.mat, '15.23')>0);
pattern = ['.*', wh, '.*flagged.*'];
files = f.mat(cellRegexpi(f.mat, pattern)>0);
erp_gen = [];isGood_gen = []; RT = [];

for i = 1:length (files)
    fprintf('Loading files: %s\n %s\n %s\n', files{i}, [wh '_' num2str(i) '_intp_csd.mat']);
    load (fullfile(fold, cell2mat (files (i))));
    load(fullfile(ERPfolder, [wh '_' num2str(i) '_intp_csd.mat']));

ImportantChans = [ 4 5 19 20 32];
isBadVolt = sq(any(maxArt(:,:,:) > 200,1)); 
isBadVoltAll = sq(any(maxArt(ImportantChans,:,:) > 100,1)); 


    isGood (isBadVolt == 1) = 0;
    isGood (isBadVoltAll == 1) = 0;
    isGood = isGood(1:864);
    RT = [RT; RS];

    %cleaning ERP
isGood_gen = [isGood_gen;isGood ];
erp_gen = nancat (3,erp_gen, erp);
end 

length (find (isGood_gen==1))
folder = 'D:\cue_task\analysis\Data\General';
fprintf('savig...');
save (fullfile(folder, [wh '_whole']),'-v7.3',...
        'isGood_gen', 'erp_gen', 'RT');
fprintf('end');
end

%eeg_data_single = single(erp_gen);

% folder = 'C:\Users\monakhvd\OneDrive - Trinity College Dublin\Data\General';
% fprintf('savig...');
% save (fullfile(folder, [wh '_whole_single']),'-v7.3',...
%         'isGood_gen', 'eeg_data_single', 'RT');
% fprintf('end');