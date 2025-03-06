clc
clear all
folder = 'G:\cue_task\analysis\Data\General';
folder1 = 'G:\cue_task\analysis\Data\cpp';

mu_files = dir(fullfile(folder, '*mu_beta.mat'));
beh_files =  dir(fullfile(folder, '*_sorted_t.mat'));
cpp_files =  dir(fullfile(folder1, '*_cpps.mat'));


beta_resp_gen = []; beta_stim_gen = []; 
cue_and_cert_gen = []; cue_and_e_gen = []; cue_m_gen = []; RT_med_gen = [];
e_or_c_gen = []; many_cond_gen = []; blocks_gen = []; cert_gen = []; pp_gen = [];lr_gen = [];

for i = 1 : length(mu_files)
load (fullfile(folder,  (mu_files(i).name)))
load (fullfile(folder,  (beh_files(i).name)))
mu_files(i).name

beta_resp = squeeze(betas.resp);
beta_stim = squeeze (betas.stim);

beta_resp_gen = nancat (2, beta_resp_gen, beta_resp); 
beta_stim_gen = nancat (2, beta_stim_gen, beta_stim); %[pp, steps, trial, contra and ipsi]
    cue_and_cert_gen = [cue_and_cert_gen; cue_and_cert]; 
    cue_and_e_gen = [cue_and_e_gen; cue_and_e]; cue_m_gen = [cue_m_gen; cue_m]; 
    e_or_c_gen = [e_or_c_gen; e_or_c]; cert_gen = [cert_gen; cert]; many_cond_gen = [many_cond_gen; many_cond];
    blocks_gen = [blocks_gen; blocks]; pp_gen = [pp_gen; pp]; lr_gen = [lr_gen; lr_resp];

load (fullfile(folder1,  (cpp_files(i).name)),'RT')

rt_med = zeros (size (RT));
 a = nanmedian (RT);
rt_med (find (RT> a)) = 1;
RT_med_gen = [RT_med_gen; rt_med];

end

c=1;
cert_and_e_gen = [];
cert_and_e_gen = NaN(size (cert_gen));

for i = 1:3
    for j = 0:1
            a = find ( cert_gen == i & e_or_c_gen == j);
            cert_and_e_gen(a) = c;
            c = c+1;
  
    end
end

fast_cue = [];
fast_cue = NaN(size (cert_gen));
c=1;
for i = 0:1
    for j = -1:1
            a = find ( RT_med_gen == i & cue_m_gen == j);
            fast_cue(a) = c;
            c = c+1;
  
    end
end

fast_cue_acc = [];
fast_cue_acc = NaN(size (cert_gen));
c=1;
for i = 0:1
    for j = 1:6
            a = find ( RT_med_gen == i & cue_and_e_gen == j);
            fast_cue_acc(a) = c;
            c = c+1;
  
    end
end

fast_cert = [];
fast_cert = NaN(size (cert_gen));
c=1;
for i = 0:1
    for j = 1:3
            a = find ( RT_med_gen == i & cert_and_e_gen == j);
            fast_cert(a) = c;
            c = c+1;
  
    end
end

fast_cert_acc = [];
fast_cert_acc = NaN(size (cert_gen));
c=1;
for i = 0:1
    for j = 1:6
            a = find ( RT_med_gen == i & cert_and_e_gen == j);
            fast_cert_acc(a) = c;
            c = c+1;
  
    end
end
%% making all matrices

beta_resp_ipsi =beta_resp_gen(:,:,1); beta_resp_contra = beta_resp_gen(:,:,2); beta_resp_diff = beta_resp_contra - beta_resp_ipsi; 
beta_stim_ipsi =beta_stim_gen(:,:,1); beta_stim_contra = beta_stim_gen(:,:,2); beta_stim_diff = beta_stim_contra - beta_stim_ipsi; 

beta_stim_diff_another = beta_stim_ipsi - beta_stim_contra;
%making it for cued hand
beta_stim_diff(:, find (cue_and_e_gen == 2 | cue_and_e_gen == 5)) =  beta_stim_diff_another(:, find (cue_and_e_gen == 2 |cue_and_e_gen == 5)) ;

pp1 = cellfun(@(x) str2num(x(2:3)), pp_gen);
cue_and_e_gen(cue_and_e_gen==0) = NaN;

dimSize = 30; 
manyCondByPP= processAndExpand(many_cond_gen, pp1, dimSize); %[TIMEPOINT, TRIAL, PP]
CertByPP = processAndExpand(cert_gen, pp1, dimSize);
CertAccByPP = processAndExpand(cert_and_e_gen, pp1, dimSize);
CertCueByPP = processAndExpand(cue_and_cert_gen, pp1, dimSize);
CueByPP = processAndExpand(cue_m_gen, pp1, dimSize);
CueAccByPP = processAndExpand(cue_and_e_gen, pp1, dimSize);
AccByPP = processAndExpand(e_or_c_gen, pp1, dimSize);
LRByPP = processAndExpand(lr_gen, pp1, dimSize);
RTMedByPP= processAndExpand(RT_med_gen, pp1, dimSize);
FastCueByPP = processAndExpand(fast_cue, pp1, dimSize);
FastCueAccByPP = processAndExpand(fast_cue_acc, pp1, dimSize);
FastCertAccByPP = processAndExpand(fast_cert_acc, pp1, dimSize);

dimSize = 65; 
StmanyCondByPP= processAndExpand(many_cond_gen, pp1, dimSize); %[TIMEPOINT, TRIAL, PP]
StCertByPP = processAndExpand(cert_gen, pp1, dimSize);
StCertAccByPP = processAndExpand(cert_and_e_gen, pp1, dimSize);
StCertCueByPP = processAndExpand(cue_and_cert_gen, pp1, dimSize);
StCueByPP = processAndExpand(cue_m_gen, pp1, dimSize);
StCueAccByPP = processAndExpand(cue_and_e_gen, pp1, dimSize);
StAccByPP = processAndExpand(e_or_c_gen, pp1, dimSize);
StLRByPP = processAndExpand(lr_gen, pp1, dimSize);
StRTMedByPP= processAndExpand(RT_med_gen, pp1, dimSize);
StFastCueByPP = processAndExpand(fast_cue, pp1, dimSize);
StFastCueAccByPP = processAndExpand(fast_cue_acc, pp1, dimSize);
StFastCertAccByPP = processAndExpand(fast_cert_acc, pp1, dimSize);

%making for resp locked

BetaRespDiffByPP = groupMeans(beta_resp_diff, 2, pp1, 'dim');
BetaRespDiffByPP = permute (BetaRespDiffByPP, [3 1 2]);


BetaStimDiffByPP = groupMeans(beta_stim_diff, 2, pp1, 'dim');
BetaStimDiffByPP = permute (BetaStimDiffByPP, [3 1 2]);


%% ipsi against contra in general

BetaStimIpsiByPP = groupMeans(beta_stim_ipsi, 2, pp1, 'dim'); %[TRIAL, PP, TIMEPOINT]
BetaStimIpsiByPP = permute (BetaStimIpsiByPP, [3 1 2]);

BetaStimContraByPP = groupMeans(beta_stim_contra, 2, pp1, 'dim');
BetaStimContraByPP = permute (BetaStimContraByPP, [3 1 2]);

BetaRespIpsiByPP = groupMeans(beta_resp_ipsi, 2, pp1, 'dim'); %[TRIAL, PP, TIMEPOINT]
BetaRespIpsiByPP = permute (BetaRespIpsiByPP, [3 1 2]);

BetaRespContraByPP = groupMeans(beta_resp_contra, 2, pp1, 'dim');
BetaRespContraByPP = permute (BetaRespContraByPP, [3 1 2]);

%
a = nanmean(permute (BetaStimIpsiByPP, [3 1 2]), 3);
b = nanmean(permute (BetaStimContraByPP, [3 1 2]), 3);

figure();

h = errorBarPlot(sq(a), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows, 'color', [.85, .325, .098]);
hold on
e = errorBarPlot(sq(b), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows, 'color', [.466, .674, .188]);
%e{1}.Color =  [.466, .674, .188];

xline(0)
legend('', 'ipsi', '', 'contra', 'Location', 'best');

% Set the title of the plot to describe its content
title(' ipsi and contra lateral');


% %% for left/right + resp locked (ipsi vs contra)
% 
% BetaRespIpsiByPPByLR = squeeze(nanmean(groupMeans(BetaRespIpsiByPP, 2, LRByPP, 'dim'), 4));
% BetaRespContraByPPByLR = squeeze(nanmean(groupMeans(BetaRespContraByPP, 2, LRByPP, 'dim'), 4));
% 
% colors_1 = [.85, .325, .098;   % Color 1
%             .466, .674, .188]; % Color 2
% 
% figure();
% colors = colors_1; set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% 
% h = errorBarPlot(sq(permute (BetaRespIpsiByPPByLR, [3 1 2])), 'area', 1, 'alpha', 0.1, 'xaxisvalues', respWindows);
% hold on
% 
% e = errorBarPlot(sq(permute (BetaRespContraByPPByLR, [3 1 2])), 'area', 1, 'alpha', 0.1, 'xaxisvalues', respWindows);
% 
% e{1}.LineStyle = "--";
% e{2}.LineStyle = "--";
% 
% xline(0)
% % Add a legend to the plot
% % Labels correspond to different error and correct conditions, with their locations specified
% legend('', 'left ipsi', '', 'right ipsi', '', 'left contra', '', 'right contra', 'Location', 'best');
% 
% % Set the title of the plot to describe its content
% title('hand for ipsi and contra lateral (response locked)');
% 
% %% for left/right + stim locked (ipsi vs contra)
% 
% BetaStimIpsipByPPByLR = squeeze(nanmean(groupMeans(BetaStimIpsiByPP, 2, StLRByPP, 'dim'), 4));
% BetaStimContrapByPPByLR = squeeze(nanmean(groupMeans(BetaStimContraByPP, 2, StLRByPP, 'dim'), 4));
% 
% colors_1 = [.85, .325, .098;   % Color 1
%             .466, .674, .188]; % Color 2
% 
% figure();
% colors = colors_1; set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% 
% h = errorBarPlot(sq(permute (BetaStimIpsipByPPByLR, [3 1 2])), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
% hold on
% 
% e = errorBarPlot(sq(permute (BetaStimContrapByPPByLR, [3 1 2])), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
% 
% e{1}.LineStyle = "--";
% e{2}.LineStyle = "--";
% 
% xline(0)
% % Add a legend to the plot
% % Labels correspond to different error and correct conditions, with their locations specified
% legend('', 'left ipsi', '', 'right ipsi', '', 'left contra', '', 'right contra', 'Location', 'best');
% 
% % Set the title of the plot to describe its content
% title('hand for ipsi and contra lateral');
% 
% 
% %% for accuracy+resp + diff contra-ipsi
% 
% BetaResDiffByPPByAcc = squeeze(nanmean(groupMeans(BetaRespDiffByPP, 2, AccByPP, 'dim'), 4));
% 
% 
% colors_1 = [.85, .325, .098;   % Color 1
%             .466, .674, .188]; % Color 2
% 
% figure();
% 
% % Define the colors for the plot using colors_5 matrix
% colors = colors_1; 
% 
% % Set the color order for the plot and specify that new plots replace existing ones in the same axes
% set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% 
% % Plot the error bars using the data cppByPPByCueE
% % The data is permuted into [trials, participants, cues] for plotting
% % 'area' specifies the plotting style, 'alpha' sets the transparency, and 'xaxisvalues' provides x-axis labels
% h = errorBarPlot(movmean(sq(permute (BetaResDiffByPPByAcc, [3 1 2])), 5, 2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', respWindows);
% 
% xline(0)
% yline(0)
% ylim ([-0.4 0.2])
% % Add a legend to the plot
% % Labels correspond to different error and correct conditions, with their locations specified
% legend('', 'error', '', 'correct',  'Location', 'best');
% 
% % Set the title of the plot to describe its content
% title('error/correct (resp locked)');
% %% for error and correct + stim locked
% 
% BetaStimDiffByPPByAcc = squeeze(nanmean(groupMeans(BetaStimDiffByPP, 2, StAccByPP, 'dim'), 4));
% 
% 
% colors_1 = [.85, .325, .098;   % Color 1
%             .466, .674, .188]; % Color 2
% 
% figure();
% 
% % Define the colors for the plot using colors_5 matrix
% colors = colors_1; 
% 
% % Set the color order for the plot and specify that new plots replace existing ones in the same axes
% set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% 
% % Plot the error bars using the data cppByPPByCueE
% 
% % The data is permuted into [trials, participants, cues] for plotting
% % 'area' specifies the plotting style, 'alpha' sets the transparency, and 'xaxisvalues' provides x-axis labels
% h = errorBarPlot(movmean(sq(permute (BetaStimDiffByPPByAcc, [3 1 2])), 5, 2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
% 
% xline(0)
% yline(0)
% ylim ([-0.3 0.3])
% x = [-1200, -900, -900, -1200]; % Define the x-coordinates
% y = get(gca, 'YLim'); % Get current y-axis limits
% patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% 
% legend('', 'error', '', 'correct', 'Location', 'best');
% 
% % Set the title of the plot to describe its content
% title('error/correct (stim locked)');
% 
% %% for cues for diff + stim locked
% 
% BetaRespDiffByPPByCue = squeeze(nanmean(groupMeans(BetaRespDiffByPP, 2, CueByPP, 'dim'), 4));
% 
% colors_1 = [.85, .325, .098;   % Color 1
%             0, .447, .741;     % Color 2
%             .466, .674, .188]; % Color 2
% 
% figure();
% 
% % Define the colors for the plot using colors_5 matrix
% colors = colors_1; 
% 
% % Set the color order for the plot and specify that new plots replace existing ones in the same axes
% set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% 
% % Plot the error bars using the data cppByPPByCueE
% % The data is permuted into [trials, participants, cues] for plotting
% % 'area' specifies the plotting style, 'alpha' sets the transparency, and 'xaxisvalues' provides x-axis labels
% h = errorBarPlot(movmean(sq(permute (BetaRespDiffByPPByCue, [3 1 2])), 5,2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', respWindows);
% 
% xline(0)
% yline(0)
% 
% % Add a legend to the plot
% % Labels correspond to different error and correct conditions, with their locations specified
% legend('', 'invalid ', '', 'neutral ','', 'valid ',   'Location', 'best');
% 
% % Set the title of the plot to describe its content
% title('cue (resp locked)');
% 
%% for cues for diff + stim locked

BetaStimDiffByPPByCue = squeeze(nanmean(groupMeans(BetaStimDiffByPP, 2, StCueByPP, 'dim'), 4));

colors_1 = [.85, .325, .098;   % Color 1
            0, .447, .741;     % Color 2
            .466, .674, .188]; % Color 2

figure();

% Define the colors for the plot using colors_5 matrix
colors = colors_1; 

% Set the color order for the plot and specify that new plots replace existing ones in the same axes
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');

% Plot the error bars using the data cppByPPByCueE
% The data is permuted into [trials, participants, cues] for plotting
% 'area' specifies the plotting style, 'alpha' sets the transparency, and 'xaxisvalues' provides x-axis labels
h = errorBarPlot(movmean(sq(permute (BetaStimDiffByPPByCue, [3 1 2])), 5,2), 'area', 1, 'alpha', 0, 'xaxisvalues', stimWindows);

xline(0)
yline(0)
x = [-1200, -900, -900, -1200]; % Define the x-coordinates
y = [-0.5 0.5]; % Get current y-axis limits
patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');


ylim ([-0.5 0.5])
xlim ([-1400 100])


  ax = gca;
    ax.Units = 'centimeters';
    ax.Position(3) = (diff(xlim) / 200); % Scale width based on xlim range

% Add a legend to the plot
% Labels correspond to different error and correct conditions, with their locations specified
legend('', 'invalid ', '', 'neutral ','', 'valid ',   'Location', 'best');

% Set the title of the plot to describe its content
title('cue for diff');
% 
% 
% %% for cues + acc for diff + resp locked
% 
% BetaRespDiffByPPByCueAcc = squeeze(nanmean(groupMeans(BetaRespDiffByPP, 2, CueAccByPP, 'dim'), 4));
% 
% colors_1 = [.85, .325, .098;  .85, .325, .098;  % Color 1
%             0, .447, .741;    0, .447, .741;   % Color 2
%            .466, .674, .188; .466, .674, .188]; % Color 2
% 
% figure();
% 
% % Define the colors for the plot using colors_5 matrix
% colors = colors_1; 
% 
% % Set the color order for the plot and specify that new plots replace existing ones in the same axes
% set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% 
% % Plot the error bars using the data cppByPPByCueE
% % The data is permuted into [trials, participants, cues] for plotting
% % 'area' specifies the plotting style, 'alpha' sets the transparency, and 'xaxisvalues' provides x-axis labels
% h = errorBarPlot(movmean(sq(permute (BetaRespDiffByPPByCueAcc, [3 1 2])), 5,2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', respWindows);
% h{1}.LineStyle = '--';h{3}.LineStyle = '--';h{5}.LineStyle = '--';
% xline(0)
% yline(0)
% ylim ([-0.8 0.8])
% xlim ([-1400 50])
% hold on;
% dummySolid = plot(nan, nan, 'k-'); % Dummy solid line for "Correct"
% dummyDashed = plot(nan, nan, 'k--'); % Dummy dashed line for "Error"
% 
% % Create custom legend entries for three main conditions and line styles
% legend([h{2}, h{4}, h{6}, dummySolid, dummyDashed], ...
%     {'Invalid', 'Neutral', 'Valid', 'Correct (Solid)', 'Error (Dashed)'}, ...
%     'Location', 'best', 'EdgeColor', 'none');
% 
% 
% title('cue (resp locked)');

%% for cues + acc for diff + stim locked

BetaStimDiffByPPByCueAcc = squeeze(nanmean(groupMeans(BetaStimDiffByPP, 2, StCueAccByPP, 'dim'), 4));

colors_1 = [.85, .325, .098;  .85, .325, .098;  % Color 1
            0, .447, .741;    0, .447, .741;   % Color 2
           .466, .674, .188; .466, .674, .188]; % Color 2

figure();

% Define the colors for the plot using colors_5 matrix
colors = colors_1; 

% Set the color order for the plot and specify that new plots replace existing ones in the same axes
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');

% Plot the error bars using the data cppByPPByCueE
% The data is permuted into [trials, participants, cues] for plotting
% 'area' specifies the plotting style, 'alpha' sets the transparency, and 'xaxisvalues' provides x-axis labels
h = errorBarPlot(movmean(sq(permute (BetaStimDiffByPPByCueAcc, [3 1 2])), 5,2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
h{1}.LineStyle = '--';h{3}.LineStyle = '--';h{5}.LineStyle = '--';
xline(0)
yline(0)
ylim ([-0.8 0.8])
xlim ([-1400 100])
x = [-1200, -900, -900, -1200]; % Define the x-coordinates
y = get(gca, 'YLim'); % Get current y-axis limits
patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

  ax = gca;
    ax.Units = 'centimeters';
    ax.Position(3) = (diff(xlim) / 200); % Scale width based on xlim range

% Create dummy lines for the legend to describe line styles
hold on;
dummySolid = plot(nan, nan, 'k-'); % Dummy solid line for "Correct"
dummyDashed = plot(nan, nan, 'k--'); % Dummy dashed line for "Error"

% Create custom legend entries for three main conditions and line styles
legend([h{2}, h{4}, h{6}, dummySolid, dummyDashed], ...
    {'Invalid', 'Neutral', 'Valid', 'Correct (Solid)', 'Error (Dashed)'}, ...
    'Location', 'best', 'EdgeColor', 'none');
% Set the title of the plot to describe its content
title('cue (stim locked)');



% %% for cues for diff + fast slow + stim locked
% 
% BetaRespDiffByPPByCueFast = squeeze(nanmean(groupMeans(BetaRespDiffByPP, 2, FastCueByPP, 'dim'), 4));
% 
% colors_1 = [.85, .325, .098;   % Color 1
%             0, .447, .741;     % Color 2
%             .466, .674, .188]; % Color 2
% 
% figure();
% 
% % Define the colors for the plot using colors_5 matrix
% colors = colors_1; 
% 
% % Set the color order for the plot and specify that new plots replace existing ones in the same axes
% set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% 
% % Plot the error bars using the data cppByPPByCueE
% % The data is permuted into [trials, participants, cues] for plotting
% % 'area' specifies the plotting style, 'alpha' sets the transparency, and 'xaxisvalues' provides x-axis labels
% h = errorBarPlot(movmean(sq(permute (BetaRespDiffByPPByCueFast, [3 1 2])), 5,2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', respWindows);
% 
% xline(0)
% yline(0)
% 
% h{4}.LineStyle = '--';h{5}.LineStyle = '--';h{6}.LineStyle = '--';
% 
% legend('', 'invalid fast', '', 'neutral fast','', 'valid fast', '', 'invalid slow', '', 'neutral slow','', 'valid slow',   'Location', 'best');
% 
% % Set the title of the plot to describe its content
% title('cue + fast/slow (resp locked)');
% 
% %% for cues for diff + fast slow + stim locked
% 
% BetaStimDiffByPPByCueFast = squeeze(nanmean(groupMeans(BetaStimDiffByPP, 2, StFastCueByPP, 'dim'), 4));
% 
% colors_1 = [.85, .325, .098;   % Color 1
%             0, .447, .741;     % Color 2
%             .466, .674, .188]; % Color 2
% 
% figure();
% 
% % Define the colors for the plot using colors_5 matrix
% colors = colors_1; 
% 
% % Set the color order for the plot and specify that new plots replace existing ones in the same axes
% set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% 
% % Plot the error bars using the data cppByPPByCueE
% % The data is permuted into [trials, participants, cues] for plotting
% % 'area' specifies the plotting style, 'alpha' sets the transparency, and 'xaxisvalues' provides x-axis labels
% h = errorBarPlot(movmean(sq(permute (BetaStimDiffByPPByCueFast, [3 1 2])), 5,2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
% 
% xline(0)
% yline(0)
% 
% h{4}.LineStyle = '--';h{5}.LineStyle = '--';h{6}.LineStyle = '--';
% 
% x = [-1200, -900, -900, -1200]; % Define the x-coordinates
% y = get(gca, 'YLim'); % Get current y-axis limits
% patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% % Add a legend to the plot
% % Labels correspond to different error and correct conditions, with their locations specified
% legend('', 'invalid fast', '', 'neutral fast','', 'valid fast', '', 'invalid slow', '', 'neutral slow','', 'valid slow',   'Location', 'best');
% 
% % Set the title of the plot to describe its content
% title('cue for diff');
% 
% % %% for cues for diff + fast  + accuracy + stim locked
% % 
% % BetaStimDiffByPPByCueFastAcc = squeeze(nanmean(groupMeans(BetaStimDiffByPP, 2, StFastCueAccByPP, 'dim'), 4));
% % 
% % colors_1 = [.85, .325, .098;   % Color 1
% %             0, .447, .741;     % Color 2
% %             .466, .674, .188;
% %              0.425, 0.1625, 0.049; 
% %              0, 0.2235, 0.3705;
% %              0.233, 0.337, 0.0944]; % Color 3
% % 
% % figure();
% % 
% % % Define the colors for the plot using colors_5 matrix
% % colors = colors_1; 
% % 
% % % Set the color order for the plot and specify that new plots replace existing ones in the same axes
% % set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% % 
% % % Plot the error bars using the data cppByPPByCueE
% % % The data is permuted into [trials, participants, cues] for plotting
% % % 'area' specifies the plotting style, 'alpha' sets the transparency, and 'xaxisvalues' provides x-axis labels
% % h = errorBarPlot(movmean(sq(permute (BetaStimDiffByPPByCueFastAcc(:, 1:6,:), [3 1 2])), 5,2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
% % 
% % xline(0)
% % yline(0)
% % h{1}.LineStyle = '--';h{2}.LineStyle = '--';h{3}.LineStyle = '--';
% % %h{7}.LineStyle = '--';h{8}.LineStyle = '--';h{9}.LineStyle = '--';
% % x = [-1200, -900, -900, -1200]; % Define the x-coordinates
% % y = get(gca, 'YLim'); % Get current y-axis limits
% % patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% % % Add a legend to the plot
% % % Labels correspond to different error and correct conditions, with their locations specified
% % legend('', 'invalid fast error', '', 'neutral fast error','', 'valid fast error',...
% %     '', 'invalid fast correct', '', 'neutral fast correct','', 'valid fast correct')%,...
% %     % '', 'invalid fast slow', '', 'neutral fast slow','', 'valid fast slow',...
% %     % '', 'invalid fast slow', '', 'neutral fast slow','', 'valid fast slow',  'Location', 'best');
% % 
% % % Set the title of the plot to describe its content
% % title('cue for diff');
% % 
% % %% for cues for diff + slow  + accuracy + stim locked
% % 
% % colors_1 = [.85, .325, .098;   % Color 1
% %             0, .447, .741;     % Color 2
% %             .466, .674, .188;
% %              0.425, 0.1625, 0.049; 
% %              0, 0.2235, 0.3705;
% %              0.233, 0.337, 0.0944]; % Color 3
% % 
% % figure();
% % 
% % % Define the colors for the plot using colors_5 matrix
% % colors = colors_1; 
% % 
% % % Set the color order for the plot and specify that new plots replace existing ones in the same axes
% % set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% % 
% % % Plot the error bars using the data cppByPPByCueE
% % % The data is permuted into [trials, participants, cues] for plotting
% % % 'area' specifies the plotting style, 'alpha' sets the transparency, and 'xaxisvalues' provides x-axis labels
% % h = errorBarPlot(movmean(sq(permute (BetaStimDiffByPPByCueFastAcc(:, 7:end,:), [3 1 2])), 5,2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
% % 
% % xline(0)
% % yline(0)
% % h{1}.LineStyle = '--';h{2}.LineStyle = '--';h{3}.LineStyle = '--';
% % 
% % x = [-1200, -900, -900, -1200]; % Define the x-coordinates
% % y = get(gca, 'YLim'); % Get current y-axis limits
% % patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% % % Add a legend to the plot
% % % Labels correspond to different error and correct conditions, with their locations specified
% % legend('', 'invalid error slow', '', 'neutral error slow','', 'valid error slow',...
% %      '', 'invalid correct slow', '', 'neutral correct slow','', 'valid correct slow',  'Location', 'best');
% % 
% % % Set the title of the plot to describe its content
% % title('cue for diff');
% 
% 
% %% for cert + resp locked
% 
% BetaRespDiffByPPByCert = squeeze(nanmean(groupMeans(BetaRespDiffByPP, 2, CertByPP, 'dim'), 4));
% 
% colors_1 = [ 0.4974,0.7088,0.9559;
%      0.5098, 0.3608, 0.7882;
%      0.7804, 0.1412, 0.4000]
% figure();
% 
% % Define the colors for the plot using colors_5 matrix
% colors = colors_1; 
% 
% % Set the color order for the plot and specify that new plots replace existing ones in the same axes
% set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% 
% h = errorBarPlot(movmean(sq(permute (BetaRespDiffByPPByCert(:,:,:), [3 1 2])), 5, 2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', respWindows);
% 
% xline(0)
% yline(0)
% 
% legend('', 'maybe ', '', 'probably  ', '', 'certain ',...
%    'Location', 'best');
% 
% % Set the title of the plot to describe its content
% title('certainty (resp locked)');
% 
% 
% %% for cert + stim locked
% 
% BetaStimDiffByPPByCert = squeeze(nanmean(groupMeans(BetaStimDiffByPP, 2, StCertByPP, 'dim'), 4));
% 
% colors_1 = [ 0.4974,0.7088,0.9559;
%      0.5098, 0.3608, 0.7882;
%      0.7804, 0.1412, 0.4000]
% figure();
% 
% % Define the colors for the plot using colors_5 matrix
% colors = colors_1; 
% 
% % Set the color order for the plot and specify that new plots replace existing ones in the same axes
% set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% 
% h = errorBarPlot(movmean(sq(permute (BetaStimDiffByPPByCert(:,:,:), [3 1 2])), 5, 2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
% 
% xline(0)
% yline(0)
% 
% x = [-1200, -900, -900, -1200]; % Define the x-coordinates
% y = get(gca, 'YLim'); % Get current y-axis limits
% patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% 
% legend('', 'maybe ', '', 'probably  ', '', 'certain ',...
%    'Location', 'best');
% 
% % Set the title of the plot to describe its content
% title('certainty (stim locked)');
% 
% 
% 
% 
%% for cert + accuracy + stim locked

BetaRespDiffByPPByCert = squeeze(nanmean(groupMeans(BetaStimDiffByPP, 2, StCertAccByPP, 'dim'), 4));

colors_1 = [ 0.4974,0.7088,0.9559; 0.4974,0.7088,0.9559;
     0.5098, 0.3608, 0.7882;0.5098, 0.3608, 0.7882;
     0.7804, 0.1412, 0.4000; 0.7804, 0.1412, 0.4000];  % Color 6
figure();

% Define the colors for the plot using colors_5 matrix
colors = colors_1; 

% Set the color order for the plot and specify that new plots replace existing ones in the same axes
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');

h = errorBarPlot(movmean(sq(permute (BetaRespDiffByPPByCert(:,1:6,:), [3 1 2])), 5, 2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
hold on

xline(0)
yline(0)
ylim ([-0.8 0.8])
xlim ([-1400 1800])
h{1}.LineStyle = "--";
h{3}.LineStyle = "--";
h{5}.LineStyle = "--";

x = [-1200, -900, -900, -1200]; % Define the x-coordinates
y = get(gca, 'YLim'); % Get current y-axis limits
patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

hold on;
dummySolid = plot(nan, nan, 'k-'); % Dummy solid line for "Correct"
dummyDashed = plot(nan, nan, 'k--'); % Dummy dashed line for "Error"

% Create custom legend entries for three main conditions and line styles
legend([h{2}, h{4}, h{6}, dummySolid, dummyDashed], ...
    {'maybe', 'probably', 'certain', 'Correct (Solid)', 'Error (Dashed)'}, ...
    'Location', 'best', 'EdgeColor', 'none');
% Set the title of the plot to describe its content
title('certainty for diff');
% 
% %% for cert + accuracy + resp locked
% 
% BetaRespDiffByPPByCert = squeeze(nanmean(groupMeans(BetaRespDiffByPP, 2, CertAccByPP, 'dim'), 4));
% 
% colors_1 = [ 0.4974,0.7088,0.9559; 0.4974,0.7088,0.9559;
%      0.5098, 0.3608, 0.7882;0.5098, 0.3608, 0.7882;
%      0.7804, 0.1412, 0.4000; 0.7804, 0.1412, 0.4000];  % Color 6
% figure();
% 
% % Define the colors for the plot using colors_5 matrix
% colors = colors_1; 
% 
% % Set the color order for the plot and specify that new plots replace existing ones in the same axes
% set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% 
% h = errorBarPlot(movmean(sq(permute (BetaRespDiffByPPByCert, [3 1 2])), 5, 2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', respWindows);
% hold on
% 
% xline(0)
% ylim ([-0.8 0.8])
% xlim ([-1400 50])
% yline(0)
% h{1}.LineStyle = "--";
% h{3}.LineStyle = "--";
% h{5}.LineStyle = "--";
% 
% 
% legend('', 'maybe error ', '', 'maybe correct ', '', 'probably error ', '', 'probably correct ', '', 'certain error ', '', 'certain correct ',...
%    'Location', 'best');
% 
% % Set the title of the plot to describe its content
% title('certainty for diff (resp)');
% 
% %% for 18 variants
% 
% BetaDiffByPPByCertCueE = squeeze(nanmean(groupMeans(BetaRespDiffByPP, 2, manyCondByPP, 'dim'), 4));
% %change order of dimentions
% BetaDiffByPPByCertCueE = movmean(sq(permute (BetaDiffByPPByCertCueE,[3,1,2])), 3,2);
% indices = 1:22;
% 
% % Define a set of colors to use for different plots
% colors_3 = [0.925, 0.6625, 0.549;  % Color 1
%             0.85, 0.325, 0.098;    % Color 2
%             0.425, 0.1625, 0.049;  % Color 3
%             0.5, 0.724, 0.8705;    % Color 4
%             0, 0.447, 0.741;       % Color 5
%             0, 0.2235, 0.3705;     % Color 6
%             0.733, 0.837, 0.594;   % Color 7
%             0.466, 0.674, 0.188;   % Color 8
%             0.233, 0.337, 0.0944]; % Color 9
% 
% % Create a figure and set up a tiled layout for the subplots
% figure();
% t = tiledlayout('flow');  % 'flow' arranges plots in a row or column based on the available space
% 
% % Create a mask identifying NaN values in the dataset
% nanMask = isnan(BetaDiffByPPByCertCueE(:, 1, :));  % Find NaN values in the first column for each participant and condition
% [row, col, page] = find(nanMask);  % Find the row and column indices where NaNs occur
% pp_to_del = [];  % Initialize an empty array to store participants to delete
% 
% % Loop over each set of 3 conditions (1:3, 4:6, etc.) to identify participants with NaNs
% for i = 1:3:18
%     % Find unique rows (participants) with NaNs in the current set of conditions
%     a = unique(row(find(col == i | col == i+1 | col == i+2)));
% 
%     % Pad with NaNs to ensure the array has 8 elements
%     a = [a; nan(15 - length(a), 1)];
% 
%     % Append the result to pp_to_del
%     pp_to_del = [pp_to_del; a'];
% end
% 
% % Rearrange the order of participants to delete, as per the specific plot order
% pp_to_del = pp_to_del([1 3 5 2 4 6], :)
% 
% % Define parameters for each plot (colors, cues, and plot titles) in a structured format
% params = {
%     {colors_3, 1:3, 'invalid+error'},          % First plot: Colors, cue indices, and plot title
%     {colors_3(4:6, :), 7:9, 'neutral+error'},  % Second plot: Colors 4-6, cue indices, and plot title
%     {colors_3(7:9, :), 13:15, 'valid+error'},  % Third plot: Colors 7-9, cue indices, and plot title
%     {colors_3, 4:6, 'invalid+correct'},        % Fourth plot: Colors, cue indices, and plot title
%     {colors_3(4:6, :), 10:12, 'neutral+correct'}, % Fifth plot
%     {colors_3(7:9, :), 16:18, 'valid+correct'} % Sixth plot
% };
% 
% % Loop through each set of parameters and create the corresponding plots
% for i = 1:length(params)
%     % Create a new tile (subplot) in the tiled layout
%     ax = nexttile(t);
% 
%     % Extract the parameters for the current plot
%     colors = params{i}{1};  % Extract the colors
%     cues = params{i}{2};    % Extract the cue indices
%     titleStr = params{i}{3};  % Extract the title for the plot
% 
%     % Identify the participants to remove based on the NaN mask (for the current plot)
%     elements_to_remove = pp_to_del(i, ~isnan(pp_to_del(i, :)));
%     modified_indices = setdiff(indices, elements_to_remove);  % Remove the participants with NaNs
% 
%     % Set the colors and plot the data
%     set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');  % Set color order for the plot
% 
%     % Plot the data using errorBarPlot with area shading and alpha transparency for each cue
%     % The data is permuted into the shape [trials, participants, cues] before plotting
%     h = errorBarPlot(BetaDiffByPPByCertCueE(modified_indices, :, cues), 'area', 1, 'alpha', 0.2, 'xaxisvalues', respWindows);
% 
%     % Add horizontal and vertical reference lines
%     yline(0);  % Horizontal line at y=0
%     xline(0);  % Vertical line at x=0
% 
%     % Set the limits for the y and x axes
%     ylim([-0.6 0.6]);  % y-axis limits from -5 to 20
%     xlim([-1500 250]);  % x-axis limits from -1500 ms to 250 ms
% 
%     % Add a legend for the plot
%     legend('', 'maybe', '', 'probably', '', 'certain', 'Location', 'best');  % Legend entries
% 
%     % Set the title for the current plot
%     title(titleStr);  % Set the title based on the parameters
% end
% 
% 
% %% for 18 variants
% 
% BetaDiffStByPPByCertCueE = squeeze(nanmean(groupMeans(BetaStimDiffByPP, 2, StmanyCondByPP, 'dim'), 4));
% %change order of dimentions
% BetaDiffStByPPByCertCueE = movmean(sq(permute (BetaDiffStByPPByCertCueE,[3,1,2])), 3,2);
% indices = 1:22;
% 
% % Define a set of colors to use for different plots
% colors_3 = [0.925, 0.6625, 0.549;  % Color 1
%             0.85, 0.325, 0.098;    % Color 2
%             0.425, 0.1625, 0.049;  % Color 3
%             0.5, 0.724, 0.8705;    % Color 4
%             0, 0.447, 0.741;       % Color 5
%             0, 0.2235, 0.3705;     % Color 6
%             0.733, 0.837, 0.594;   % Color 7
%             0.466, 0.674, 0.188;   % Color 8
%             0.233, 0.337, 0.0944]; % Color 9
% 
% % Create a figure and set up a tiled layout for the subplots
% figure();
% t = tiledlayout('flow');  % 'flow' arranges plots in a row or column based on the available space
% 
% % Create a mask identifying NaN values in the dataset
% nanMask = isnan(BetaDiffStByPPByCertCueE(:, 1, :));  % Find NaN values in the first column for each participant and condition
% [row, col, page] = find(nanMask);  % Find the row and column indices where NaNs occur
% pp_to_del = [];  % Initialize an empty array to store participants to delete
% 
% % Loop over each set of 3 conditions (1:3, 4:6, etc.) to identify participants with NaNs
% for i = 1:3:18
%     % Find unique rows (participants) with NaNs in the current set of conditions
%     a = unique(row(find(col == i | col == i+1 | col == i+2)));
% 
%     % Pad with NaNs to ensure the array has 8 elements
%     a = [a; nan(15 - length(a), 1)];
% 
%     % Append the result to pp_to_del
%     pp_to_del = [pp_to_del; a'];
% end
% 
% % Rearrange the order of participants to delete, as per the specific plot order
% pp_to_del = pp_to_del([1 3 5 2 4 6], :)
% 
% % Define parameters for each plot (colors, cues, and plot titles) in a structured format
% params = {
%     {colors_3, 1:3, 'invalid+error'},          % First plot: Colors, cue indices, and plot title
%     {colors_3(4:6, :), 7:9, 'neutral+error'},  % Second plot: Colors 4-6, cue indices, and plot title
%     {colors_3(7:9, :), 13:15, 'valid+error'},  % Third plot: Colors 7-9, cue indices, and plot title
%     {colors_3, 4:6, 'invalid+correct'},        % Fourth plot: Colors, cue indices, and plot title
%     {colors_3(4:6, :), 10:12, 'neutral+correct'}, % Fifth plot
%     {colors_3(7:9, :), 16:18, 'valid+correct'} % Sixth plot
% };
% 
% % Loop through each set of parameters and create the corresponding plots
% for i = 1:length(params)
%     % Create a new tile (subplot) in the tiled layout
%     ax = nexttile(t);
% 
%     % Extract the parameters for the current plot
%     colors = params{i}{1};  % Extract the colors
%     cues = params{i}{2};    % Extract the cue indices
%     titleStr = params{i}{3};  % Extract the title for the plot
% 
%     % Identify the participants to remove based on the NaN mask (for the current plot)
%     elements_to_remove = pp_to_del(i, ~isnan(pp_to_del(i, :)));
%     modified_indices = setdiff(indices, elements_to_remove);  % Remove the participants with NaNs
% 
%     % Set the colors and plot the data
%     set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');  % Set color order for the plot
% 
%     % Plot the data using errorBarPlot with area shading and alpha transparency for each cue
%     % The data is permuted into the shape [trials, participants, cues] before plotting
%     h = errorBarPlot(BetaDiffStByPPByCertCueE(modified_indices, :, cues), 'area', 1, 'alpha', 0.2, 'xaxisvalues', stimWindows);
% 
%     % Add horizontal and vertical reference lines
%     yline(0);  % Horizontal line at y=0
%     xline(0);  % Vertical line at x=0
% 
%     % Set the limits for the y and x axes
%     ylim([-0.6 0.6]);  % y-axis limits from -5 to 20
% 
% 
% 
%     x = [-1200, -900, -900, -1200]; % Define the x-coordinates
% y = get(gca, 'YLim'); % Get current y-axis limits
% patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% 
% 
%     % Add a legend for the plot
%     legend('', 'maybe', '', 'probably', '', 'certain', 'Location', 'best');  % Legend entries
% 
%     % Set the title for the current plot
%     title(titleStr);  % Set the title based on the parameters
% end
% 
% 
% 
% 
% 
% 
% 
% %% for 18 variants for ipsi
% 
% BetaIpsiByPPByCertCueE = squeeze(nanmean(groupMeans(BetaRespIpsiByPP, 2, manyCondByPP, 'dim'), 4));
% %change order of dimentions
% BetaIpsiByPPByCertCueE = movmean(sq(permute (BetaIpsiByPPByCertCueE,[3,1,2])), 3,2);
% indices = 1:22;
% 
% % Define a set of colors to use for different plots
% colors_3 = [0.925, 0.6625, 0.549;  % Color 1
%             0.85, 0.325, 0.098;    % Color 2
%             0.425, 0.1625, 0.049;  % Color 3
%             0.5, 0.724, 0.8705;    % Color 4
%             0, 0.447, 0.741;       % Color 5
%             0, 0.2235, 0.3705;     % Color 6
%             0.733, 0.837, 0.594;   % Color 7
%             0.466, 0.674, 0.188;   % Color 8
%             0.233, 0.337, 0.0944]; % Color 9
% 
% % Create a figure and set up a tiled layout for the subplots
% figure();
% t = tiledlayout('flow');  % 'flow' arranges plots in a row or column based on the available space
% 
% % Create a mask identifying NaN values in the dataset
% nanMask = isnan(BetaIpsiByPPByCertCueE(:, 1, :));  % Find NaN values in the first column for each participant and condition
% [row, col, page] = find(nanMask);  % Find the row and column indices where NaNs occur
% pp_to_del = [];  % Initialize an empty array to store participants to delete
% 
% % Loop over each set of 3 conditions (1:3, 4:6, etc.) to identify participants with NaNs
% for i = 1:3:18
%     % Find unique rows (participants) with NaNs in the current set of conditions
%     a = unique(row(find(col == i | col == i+1 | col == i+2)));
% 
%     % Pad with NaNs to ensure the array has 8 elements
%     a = [a; nan(15 - length(a), 1)];
% 
%     % Append the result to pp_to_del
%     pp_to_del = [pp_to_del; a'];
% end
% 
% % Rearrange the order of participants to delete, as per the specific plot order
% pp_to_del = pp_to_del([1 3 5 2 4 6], :)
% 
% % Define parameters for each plot (colors, cues, and plot titles) in a structured format
% params = {
%     {colors_3, 1:3, 'invalid+error +ipsi'},          % First plot: Colors, cue indices, and plot title
%     {colors_3(4:6, :), 7:9, 'neutral+error+ipsi'},  % Second plot: Colors 4-6, cue indices, and plot title
%     {colors_3(7:9, :), 13:15, 'valid+error+ipsi'},  % Third plot: Colors 7-9, cue indices, and plot title
%     {colors_3, 4:6, 'invalid+correct+ipsi'},        % Fourth plot: Colors, cue indices, and plot title
%     {colors_3(4:6, :), 10:12, 'neutral+correct+ipsi'}, % Fifth plot
%     {colors_3(7:9, :), 16:18, 'valid+correct+ipsi'} % Sixth plot
% };
% 
% % Loop through each set of parameters and create the corresponding plots
% for i = 1:length(params)
%     % Create a new tile (subplot) in the tiled layout
%     ax = nexttile(t);
% 
%     % Extract the parameters for the current plot
%     colors = params{i}{1};  % Extract the colors
%     cues = params{i}{2};    % Extract the cue indices
%     titleStr = params{i}{3};  % Extract the title for the plot
% 
%     % Identify the participants to remove based on the NaN mask (for the current plot)
%     elements_to_remove = pp_to_del(i, ~isnan(pp_to_del(i, :)));
%     modified_indices = setdiff(indices, elements_to_remove);  % Remove the participants with NaNs
% 
%     % Set the colors and plot the data
%     set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');  % Set color order for the plot
% 
%     % Plot the data using errorBarPlot with area shading and alpha transparency for each cue
%     % The data is permuted into the shape [trials, participants, cues] before plotting
%     h = errorBarPlot(BetaIpsiByPPByCertCueE(modified_indices, :, cues), 'area', 1, 'alpha', 0.2, 'xaxisvalues', respWindows);
% 
%     % Add horizontal and vertical reference lines
%     yline(0);  % Horizontal line at y=0
%     xline(0);  % Vertical line at x=0
%     % Set the limits for the y and x axes
%     ylim([4.3 6.7]);  % y-axis limits from -5 to 20
%     xlim([-1500 250]);  % x-axis limits from -1500 ms to 250 ms
% 
%     % Add a legend for the plot
%     legend('', 'maybe', '', 'probably', '', 'certain', 'Location', 'best');  % Legend entries
% 
%     % Set the title for the current plot
%     title(titleStr);  % Set the title based on the parameters
% end
% 
% 
% %% for 18 variants (stim) for ipsi
% 
% BetaIpsiStByPPByCertCueE = squeeze(nanmean(groupMeans(BetaStimIpsiByPP, 2, StmanyCondByPP, 'dim'), 4));
% %change order of dimentions
% BetaIpsiStByPPByCertCueE = movmean(sq(permute (BetaIpsiStByPPByCertCueE,[3,1,2])), 3,2);
% indices = 1:22;
% 
% % Define a set of colors to use for different plots
% colors_3 = [0.925, 0.6625, 0.549;  % Color 1
%             0.85, 0.325, 0.098;    % Color 2
%             0.425, 0.1625, 0.049;  % Color 3
%             0.5, 0.724, 0.8705;    % Color 4
%             0, 0.447, 0.741;       % Color 5
%             0, 0.2235, 0.3705;     % Color 6
%             0.733, 0.837, 0.594;   % Color 7
%             0.466, 0.674, 0.188;   % Color 8
%             0.233, 0.337, 0.0944]; % Color 9
% 
% % Create a figure and set up a tiled layout for the subplots
% figure();
% t = tiledlayout('flow');  % 'flow' arranges plots in a row or column based on the available space
% 
% % Create a mask identifying NaN values in the dataset
% nanMask = isnan(BetaIpsiStByPPByCertCueE(:, 1, :));  % Find NaN values in the first column for each participant and condition
% [row, col, page] = find(nanMask);  % Find the row and column indices where NaNs occur
% pp_to_del = [];  % Initialize an empty array to store participants to delete
% 
% % Loop over each set of 3 conditions (1:3, 4:6, etc.) to identify participants with NaNs
% for i = 1:3:18
%     % Find unique rows (participants) with NaNs in the current set of conditions
%     a = unique(row(find(col == i | col == i+1 | col == i+2)));
% 
%     % Pad with NaNs to ensure the array has 8 elements
%     a = [a; nan(15 - length(a), 1)];
% 
%     % Append the result to pp_to_del
%     pp_to_del = [pp_to_del; a'];
% end
% 
% % Rearrange the order of participants to delete, as per the specific plot order
% pp_to_del = pp_to_del([1 3 5 2 4 6], :)
% 
% % Define parameters for each plot (colors, cues, and plot titles) in a structured format
% params = {
%     {colors_3, 1:3, 'invalid+error+ipsi'},          % First plot: Colors, cue indices, and plot title
%     {colors_3(4:6, :), 7:9, 'neutral+error+ipsi'},  % Second plot: Colors 4-6, cue indices, and plot title
%     {colors_3(7:9, :), 13:15, 'valid+error+ipsi'},  % Third plot: Colors 7-9, cue indices, and plot title
%     {colors_3, 4:6, 'invalid+correct+ipsi'},        % Fourth plot: Colors, cue indices, and plot title
%     {colors_3(4:6, :), 10:12, 'neutral+correct+ipsi'}, % Fifth plot
%     {colors_3(7:9, :), 16:18, 'valid+correct+ipsi'} % Sixth plot
% };
% 
% % Loop through each set of parameters and create the corresponding plots
% for i = 1:length(params)
%     % Create a new tile (subplot) in the tiled layout
%     ax = nexttile(t);
% 
%     % Extract the parameters for the current plot
%     colors = params{i}{1};  % Extract the colors
%     cues = params{i}{2};    % Extract the cue indices
%     titleStr = params{i}{3};  % Extract the title for the plot
% 
%     % Identify the participants to remove based on the NaN mask (for the current plot)
%     elements_to_remove = pp_to_del(i, ~isnan(pp_to_del(i, :)));
%     modified_indices = setdiff(indices, elements_to_remove);  % Remove the participants with NaNs
% 
%     % Set the colors and plot the data
%     set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');  % Set color order for the plot
% 
%     % Plot the data using errorBarPlot with area shading and alpha transparency for each cue
%     % The data is permuted into the shape [trials, participants, cues] before plotting
%     h = errorBarPlot(BetaIpsiStByPPByCertCueE(modified_indices, :, cues), 'area', 1, 'alpha', 0.2, 'xaxisvalues', stimWindows);
% 
%     % Add horizontal and vertical reference lines
%     yline(0);  % Horizontal line at y=0
%     xline(0);  % Vertical line at x=0
% 
%     % Set the limits for the y and x axes
%     ylim([4.3 7.2]);  % y-axis limits from -5 to 20
% 
% 
% 
%     x = [-1200, -900, -900, -1200]; % Define the x-coordinates
% y = get(gca, 'YLim'); % Get current y-axis limits
% patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% 
% 
%     % Add a legend for the plot
%     legend('', 'maybe', '', 'probably', '', 'certain', 'Location', 'best');  % Legend entries
% 
%     % Set the title for the current plot
%     title(titleStr);  % Set the title based on the parameters
% end
% 
% %% for 18 variants for contra
% 
% BetaContraByPPByCertCueE = squeeze(nanmean(groupMeans(BetaRespContraByPP, 2, manyCondByPP, 'dim'), 4));
% %change order of dimentions
% BetaContraByPPByCertCueE = movmean(sq(permute (BetaContraByPPByCertCueE,[3,1,2])), 3,2);
% indices = 1:22;
% 
% % Define a set of colors to use for different plots
% colors_3 = [0.925, 0.6625, 0.549;  % Color 1
%             0.85, 0.325, 0.098;    % Color 2
%             0.425, 0.1625, 0.049;  % Color 3
%             0.5, 0.724, 0.8705;    % Color 4
%             0, 0.447, 0.741;       % Color 5
%             0, 0.2235, 0.3705;     % Color 6
%             0.733, 0.837, 0.594;   % Color 7
%             0.466, 0.674, 0.188;   % Color 8
%             0.233, 0.337, 0.0944]; % Color 9
% 
% % Create a figure and set up a tiled layout for the subplots
% figure();
% t = tiledlayout('flow');  % 'flow' arranges plots in a row or column based on the available space
% 
% % Create a mask identifying NaN values in the dataset
% nanMask = isnan(BetaContraByPPByCertCueE(:, 1, :));  % Find NaN values in the first column for each participant and condition
% [row, col, page] = find(nanMask);  % Find the row and column indices where NaNs occur
% pp_to_del = [];  % Initialize an empty array to store participants to delete
% 
% % Loop over each set of 3 conditions (1:3, 4:6, etc.) to identify participants with NaNs
% for i = 1:3:18
%     % Find unique rows (participants) with NaNs in the current set of conditions
%     a = unique(row(find(col == i | col == i+1 | col == i+2)));
% 
%     % Pad with NaNs to ensure the array has 8 elements
%     a = [a; nan(15 - length(a), 1)];
% 
%     % Append the result to pp_to_del
%     pp_to_del = [pp_to_del; a'];
% end
% 
% % Rearrange the order of participants to delete, as per the specific plot order
% pp_to_del = pp_to_del([1 3 5 2 4 6], :)
% 
% % Define parameters for each plot (colors, cues, and plot titles) in a structured format
% params = {
%     {colors_3, 1:3, 'invalid+error +contra'},          % First plot: Colors, cue indices, and plot title
%     {colors_3(4:6, :), 7:9, 'neutral+error+contra'},  % Second plot: Colors 4-6, cue indices, and plot title
%     {colors_3(7:9, :), 13:15, 'valid+error+contra'},  % Third plot: Colors 7-9, cue indices, and plot title
%     {colors_3, 4:6, 'invalid+correct+contra'},        % Fourth plot: Colors, cue indices, and plot title
%     {colors_3(4:6, :), 10:12, 'neutral+correct+contra'}, % Fifth plot
%     {colors_3(7:9, :), 16:18, 'valid+correct+contra'} % Sixth plot
% };
% 
% % Loop through each set of parameters and create the corresponding plots
% for i = 1:length(params)
%     % Create a new tile (subplot) in the tiled layout
%     ax = nexttile(t);
% 
%     % Extract the parameters for the current plot
%     colors = params{i}{1};  % Extract the colors
%     cues = params{i}{2};    % Extract the cue indices
%     titleStr = params{i}{3};  % Extract the title for the plot
% 
%     % Identify the participants to remove based on the NaN mask (for the current plot)
%     elements_to_remove = pp_to_del(i, ~isnan(pp_to_del(i, :)));
%     modified_indices = setdiff(indices, elements_to_remove);  % Remove the participants with NaNs
% 
%     % Set the colors and plot the data
%     set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');  % Set color order for the plot
% 
%     % Plot the data using errorBarPlot with area shading and alpha transparency for each cue
%     % The data is permuted into the shape [trials, participants, cues] before plotting
%     h = errorBarPlot(BetaContraByPPByCertCueE(modified_indices, :, cues), 'area', 1, 'alpha', 0.2, 'xaxisvalues', respWindows);
% 
%     % Add horizontal and vertical reference lines
%     yline(0);  % Horizontal line at y=0
%     xline(0);  % Vertical line at x=0
%     % Set the limits for the y and x axes
%     ylim([4.3 6.7]);  % y-axis limits from -5 to 20
%     xlim([-1500 250]);  % x-axis limits from -1500 ms to 250 ms
% 
%     % Add a legend for the plot
%     legend('', 'maybe', '', 'probably', '', 'certain', 'Location', 'best');  % Legend entries
% 
%     % Set the title for the current plot
%     title(titleStr);  % Set the title based on the parameters
% end
% 
% 
% %% for 18 variants (stim) +contra
% 
% BetaContraStByPPByCertCueE = squeeze(nanmean(groupMeans(BetaStimContraByPP, 2, StmanyCondByPP, 'dim'), 4));
% %change order of dimentions
% BetaContraStByPPByCertCueE = movmean(sq(permute (BetaContraStByPPByCertCueE,[3,1,2])), 3,2);
% indices = 1:22;
% 
% % Define a set of colors to use for different plots
% colors_3 = [0.925, 0.6625, 0.549;  % Color 1
%             0.85, 0.325, 0.098;    % Color 2
%             0.425, 0.1625, 0.049;  % Color 3
%             0.5, 0.724, 0.8705;    % Color 4
%             0, 0.447, 0.741;       % Color 5
%             0, 0.2235, 0.3705;     % Color 6
%             0.733, 0.837, 0.594;   % Color 7
%             0.466, 0.674, 0.188;   % Color 8
%             0.233, 0.337, 0.0944]; % Color 9
% 
% % Create a figure and set up a tiled layout for the subplots
% figure();
% t = tiledlayout('flow');  % 'flow' arranges plots in a row or column based on the available space
% 
% % Create a mask identifying NaN values in the dataset
% nanMask = isnan(BetaContraStByPPByCertCueE(:, 1, :));  % Find NaN values in the first column for each participant and condition
% [row, col, page] = find(nanMask);  % Find the row and column indices where NaNs occur
% pp_to_del = [];  % Initialize an empty array to store participants to delete
% 
% % Loop over each set of 3 conditions (1:3, 4:6, etc.) to identify participants with NaNs
% for i = 1:3:18
%     % Find unique rows (participants) with NaNs in the current set of conditions
%     a = unique(row(find(col == i | col == i+1 | col == i+2)));
% 
%     % Pad with NaNs to ensure the array has 8 elements
%     a = [a; nan(15 - length(a), 1)];
% 
%     % Append the result to pp_to_del
%     pp_to_del = [pp_to_del; a'];
% end
% 
% % Rearrange the order of participants to delete, as per the specific plot order
% pp_to_del = pp_to_del([1 3 5 2 4 6], :)
% 
% % Define parameters for each plot (colors, cues, and plot titles) in a structured format
% params = {
%     {colors_3, 1:3, 'invalid+error+contra'},          % First plot: Colors, cue indices, and plot title
%     {colors_3(4:6, :), 7:9, 'neutral+error+contra'},  % Second plot: Colors 4-6, cue indices, and plot title
%     {colors_3(7:9, :), 13:15, 'valid+error+contra'},  % Third plot: Colors 7-9, cue indices, and plot title
%     {colors_3, 4:6, 'invalid+correct+contra'},        % Fourth plot: Colors, cue indices, and plot title
%     {colors_3(4:6, :), 10:12, 'neutral+correct+contra'}, % Fifth plot
%     {colors_3(7:9, :), 16:18, 'valid+correct+contra'} % Sixth plot
% };
% 
% % Loop through each set of parameters and create the corresponding plots
% for i = 1:length(params)
%     % Create a new tile (subplot) in the tiled layout
%     ax = nexttile(t);
% 
%     % Extract the parameters for the current plot
%     colors = params{i}{1};  % Extract the colors
%     cues = params{i}{2};    % Extract the cue indices
%     titleStr = params{i}{3};  % Extract the title for the plot
% 
%     % Identify the participants to remove based on the NaN mask (for the current plot)
%     elements_to_remove = pp_to_del(i, ~isnan(pp_to_del(i, :)));
%     modified_indices = setdiff(indices, elements_to_remove);  % Remove the participants with NaNs
% 
%     % Set the colors and plot the data
%     set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');  % Set color order for the plot
% 
%     % Plot the data using errorBarPlot with area shading and alpha transparency for each cue
%     % The data is permuted into the shape [trials, participants, cues] before plotting
%     h = errorBarPlot(BetaContraStByPPByCertCueE(modified_indices, :, cues), 'area', 1, 'alpha', 0.2, 'xaxisvalues', stimWindows);
% 
%     % Add horizontal and vertical reference lines
%     yline(0);  % Horizontal line at y=0
%     xline(0);  % Vertical line at x=0
% 
%     % Set the limits for the y and x axes
%     ylim([4.3 7.2]);  % y-axis limits from -5 to 20
% 
% 
% 
%     x = [-1200, -900, -900, -1200]; % Define the x-coordinates
% y = get(gca, 'YLim'); % Get current y-axis limits
% patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% 
% 
%     % Add a legend for the plot
%     legend('', 'maybe', '', 'probably', '', 'certain', 'Location', 'best');  % Legend entries
% 
%     % Set the title for the current plot
%     title(titleStr);  % Set the title based on the parameters
% end
% 
% %% for cues + acc for ipsi + resp locked
% 
% BetaRespIpsiByPPByCueAcc = squeeze(nanmean(groupMeans(BetaRespIpsiByPP, 2, CueAccByPP, 'dim'), 4));
% 
% colors_1 = [.85, .325, .098;  .85, .325, .098;  % Color 1
%             0, .447, .741;    0, .447, .741;   % Color 2
%            .466, .674, .188; .466, .674, .188]; % Color 2
% 
% figure();
% 
% % Define the colors for the plot using colors_5 matrix
% colors = colors_1; 
% 
% % Set the color order for the plot and specify that new plots replace existing ones in the same axes
% set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% 
% % Plot the error bars using the data cppByPPByCueE
% % The data is permuted into [trials, participants, cues] for plotting
% % 'area' specifies the plotting style, 'alpha' sets the transparency, and 'xaxisvalues' provides x-axis labels
% h = errorBarPlot(movmean(sq(permute (BetaRespIpsiByPPByCueAcc, [3 1 2])), 5,2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', respWindows);
% h{1}.LineStyle = '--';h{3}.LineStyle = '--';h{5}.LineStyle = '--';
% xline(0)
% ylim ([4.6 6.6])
% % yline(0)
% 
% legend('', 'invalid error', '', 'invalid correct', '', 'neutral error','', 'neutral correct','', 'valid error','', 'valid correct',   'Location', 'best');
% 
% 
% title('cue ipsi(resp locked)');
% 
% %% for cues + acc for ipsi + resp locked
% 
% BetaRespContraByPPByCueAcc = squeeze(nanmean(groupMeans(BetaRespContraByPP, 2, CueAccByPP, 'dim'), 4));
% 
% colors_1 = [.85, .325, .098;  .85, .325, .098;  % Color 1
%             0, .447, .741;    0, .447, .741;   % Color 2
%            .466, .674, .188; .466, .674, .188]; % Color 2
% 
% figure();
% 
% % Define the colors for the plot using colors_5 matrix
% colors = colors_1; 
% 
% % Set the color order for the plot and specify that new plots replace existing ones in the same axes
% set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% 
% % Plot the error bars using the data cppByPPByCueE
% % The data is permuted into [trials, participants, cues] for plotting
% % 'area' specifies the plotting style, 'alpha' sets the transparency, and 'xaxisvalues' provides x-axis labels
% h = errorBarPlot(movmean(sq(permute (BetaRespContraByPPByCueAcc, [3 1 2])), 5,2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', respWindows);
% h{1}.LineStyle = '--';h{3}.LineStyle = '--';h{5}.LineStyle = '--';
% xline(0)
% ylim ([4.6 6.6])
% % yline(0)
% 
% legend('', 'invalid error', '', 'invalid correct', '', 'neutral error','', 'neutral correct','', 'valid error','', 'valid correct',   'Location', 'best');
% 
% 
% title('cue contra(resp locked)');
% %% for cues + acc for ipsi+ stim locked
% 
% BetaStimIpsiByPPByCueAcc = squeeze(nanmean(groupMeans(BetaStimIpsiByPP, 2, StCueAccByPP, 'dim'), 4));
% 
% colors_1 = [.85, .325, .098;  .85, .325, .098;  % Color 1
%             0, .447, .741;    0, .447, .741;   % Color 2
%            .466, .674, .188; .466, .674, .188]; % Color 2
% 
% figure();
% 
% % Define the colors for the plot using colors_5 matrix
% colors = colors_1; 
% 
% % Set the color order for the plot and specify that new plots replace existing ones in the same axes
% set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% 
% % Plot the error bars using the data cppByPPByCueE
% % The data is permuted into [trials, participants, cues] for plotting
% % 'area' specifies the plotting style, 'alpha' sets the transparency, and 'xaxisvalues' provides x-axis labels
% h = errorBarPlot(movmean(sq(permute (BetaStimIpsiByPPByCueAcc, [3 1 2])), 5,2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
% h{1}.LineStyle = '--';h{3}.LineStyle = '--';h{5}.LineStyle = '--';
% xline(0)
% ylim([4.5 7])
% %yline(0)
% x = [-1200, -900, -900, -1200]; % Define the x-coordinates
% y = get(gca, 'YLim'); % Get current y-axis limits
% patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% % Add a legend to the plot
% % Labels correspond to different error and correct conditions, with their locations specified
% legend('', 'invalid error', '', 'invalid correct', '', 'neutral error','', 'neutral correct','', 'valid error','', 'valid correct',   'Location', 'best');
% 
% % Set the title of the plot to describe its content
% title('cue ipsi (stim locked)');
% 
% 
% %% for cues + acc for ipsi+ stim locked
% 
% BetaStimCntraByPPByCueAcc = squeeze(nanmean(groupMeans(BetaStimContraByPP, 2, StCueAccByPP, 'dim'), 4));
% 
% colors_1 = [.85, .325, .098;  .85, .325, .098;  % Color 1
%             0, .447, .741;    0, .447, .741;   % Color 2
%            .466, .674, .188; .466, .674, .188]; % Color 2
% 
% figure();
% 
% % Define the colors for the plot using colors_5 matrix
% colors = colors_1; 
% 
% % Set the color order for the plot and specify that new plots replace existing ones in the same axes
% set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
% 
% % Plot the error bars using the data cppByPPByCueE
% % The data is permuted into [trials, participants, cues] for plotting
% % 'area' specifies the plotting style, 'alpha' sets the transparency, and 'xaxisvalues' provides x-axis labels
% h = errorBarPlot(movmean(sq(permute (BetaStimCntraByPPByCueAcc, [3 1 2])), 5,2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
% h{1}.LineStyle = '--';h{3}.LineStyle = '--';h{5}.LineStyle = '--';
% xline(0)
% ylim([4.5 7])
% %yline(0)
% x = [-1200, -900, -900, -1200]; % Define the x-coordinates
% y = get(gca, 'YLim'); % Get current y-axis limits
% patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% % Add a legend to the plot
% % Labels correspond to different error and correct conditions, with their locations specified
% legend('', 'invalid error', '', 'invalid correct', '', 'neutral error','', 'neutral correct','', 'valid error','', 'valid correct',   'Location', 'best');
% 
% % Set the title of the plot to describe its content
% title('cue comtra(stim locked)');