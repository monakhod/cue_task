clear all clc
folder = 'G:\cue_task\analysis\Data\General';
mat_files = dir(fullfile(folder, '*_difference_map_8to30.mat')); %getting files with mean for 8-30 Hz
chanInfo = load('chanlocsBioSemi128.mat'); % note this was made by calling>> readlocs('cap128.loc') , and, since the locations for the 8 externals (129:136) are all (1,0,0), getting rid of those by calling chanlocs = chanlocs(1:128)
eeg.chanlocs = chanInfo.chanlocs;
eeg.chanNames = {eeg.chanlocs.labels}';
 map_plot_gen = [];

cppChanNames = {'D18','D19','D28', 'B22','B21', 'B18'}';
    [~, cppChanInds] = ismember(cppChanNames, eeg.chanNames);

 for i = 1:length (mat_files)
    fprintf('Loading files: %s\n ', mat_files(i).name);
    load (fullfile(folder,  (mat_files(i).name)));

     map_plot_gen = [map_plot_gen, difference];
 end

n = {'P01', 'P02', 'P03', 'P04', 'P05', 'P06',...
    'P07', 'P08', 'P10', ...
    'P13', 'P14', 'P15', 'P16', 'P17', 'P18' 'P19' 'P20' 'P21' 'P22'};
if ~exist('topoplot','file')
    eeglab nogui;
end

 figure;

for i = 1:19
    subplot(5, 4, i);
    topoplot(map_plot_gen(:, i), eeg.chanlocs, 'electrodes','off','colormap',crameri('vik'),...
        'emarker',{'.','k',10,1},'emarker2',{col(cppChanInds), '.','y',5,1},...
        'mapLimits', [-0.6 0.6]);
    colorbar;
    title(['Difference for participant ' n{i}]);
end

figure;
a = nanmean (map_plot_gen, 2);
    topoplot(a, eeg.chanlocs, 'electrodes','on','colormap',crameri('vik'),...
    'emarker',{'.','k',5,1},'emarker2',{col(cppChanInds), '.','y',7,1});
colorbar
% % pick two clusters (left + right hemisphere)
chanClusters = reshape(cppChanNames,3,2);
[~, chansInds] = ismember(col(chanClusters), eeg.chanNames);