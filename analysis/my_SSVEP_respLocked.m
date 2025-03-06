clc
clear all
folder = 'D:\cue_task\analysis\Data\General';
folder1 = 'D:\cue_task\analysis\Data\cpp';

ssvep_files = dir(fullfile(folder, '*_SSVEP.mat'));
beh_files =  dir(fullfile(folder, '*_sorted_t.mat'));
cpp_files =  dir(fullfile(folder1, '*_cpps.mat'));

 SSVEP_tweny_gen = []; SSVEP_twenyFive_gen = [];
cue_and_cert_gen = []; cue_and_e_gen = []; cue_m_gen = []; RT_med_gen = [];
e_or_c_gen = []; many_cond_gen = []; blocks_gen = []; cert_gen = []; pp_gen = [];lr_gen = [];
freq_correct_gen = []; freq_ch_gen = []; RT_gen = [];

for i = 1 : length(ssvep_files)
load (fullfile(folder,  (ssvep_files(i).name)))
load (fullfile(folder,  (beh_files(i).name)))
ssvep_files(i).name


SSVEP_twenyresp =  squeeze (ssvep.twenty_resp);
SSVEP_twenyFiveresp = squeeze (ssvep.twentyFive_resp);

SSVEP_tweny_gen = nancat (2, SSVEP_tweny_gen, SSVEP_twenyresp); %[pp, steps, trial]
SSVEP_twenyFive_gen = nancat (2, SSVEP_twenyFive_gen, SSVEP_twenyFiveresp); %[pp, steps, trial]

    cue_and_cert_gen = [cue_and_cert_gen; cue_and_cert]; 
    cue_and_e_gen = [cue_and_e_gen; cue_and_e]; cue_m_gen = [cue_m_gen; cue_m]; 
    e_or_c_gen = [e_or_c_gen; e_or_c]; cert_gen = [cert_gen; cert]; many_cond_gen = [many_cond_gen; many_cond];
    blocks_gen = [blocks_gen; blocks]; pp_gen = [pp_gen; pp]; lr_gen = [lr_gen; lr_resp];
    freq_correct_gen = [freq_correct_gen; freq_correct]; freq_ch_gen = [freq_ch_gen; freq_ch];

load (fullfile(folder1,  (cpp_files(i).name)),'RT')
 RT_gen = [RT_gen; RT];
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


session_gen = NaN(size (cert_gen));
a = find(ismember(blocks_gen, 1:6));
session_gen(a) = 1;
a = find(ismember(blocks_gen, 7:12));
session_gen(a) = 2;
a = find(ismember(blocks_gen, 13:18));
session_gen(a) = 3;

respWindows = -1400:50:50;

%% making all matrices
SSVEP_target = SSVEP_tweny_gen;
SSVEP_target(:, freq_correct_gen==2) = SSVEP_twenyFive_gen(:, freq_correct_gen==2);

SSVEP_non_target = SSVEP_tweny_gen;
SSVEP_non_target(:, freq_correct_gen==1) = SSVEP_twenyFive_gen(:, freq_correct_gen==1);


SSVEP_diff = SSVEP_target - SSVEP_non_target;

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
RFrCorrByPP = processAndExpand(freq_correct_gen, pp1, dimSize);
RSessionByPP = processAndExpand(session_gen, pp1, dimSize);
%making for resp locked

SSVEPTwenyByPP = groupMeans(SSVEP_tweny_gen, 2, pp1, 'dim');% [trial, pp, step]
SSVEPTwenyByPP = permute (SSVEPTwenyByPP, [3 1 2]);


SSVEPTwenyFiveByPP = groupMeans(SSVEP_twenyFive_gen, 2, pp1, 'dim');% [trial, pp, step]
SSVEPTwenyFiveByPP = permute (SSVEPTwenyFiveByPP, [3 1 2]);

SSVEP_targetByPP = groupMeans(SSVEP_target, 2, pp1, 'dim');% [trial, pp, step]
SSVEP_targetByPP = permute (SSVEP_targetByPP, [3 1 2]);

SSVEP_nontargetByPP = groupMeans(SSVEP_non_target, 2, pp1, 'dim');% [trial, pp, step]
SSVEP_nontargetByPP = permute (SSVEP_nontargetByPP, [3 1 2]);


SSVEPdiffByPP = groupMeans(SSVEP_diff, 2, pp1, 'dim');% [trial, pp, step]
SSVEPdiffByPP = permute (SSVEPdiffByPP, [3 1 2]);
%% ploting plots

SSVEPTwenyByPPByHZCorr = squeeze(nanmean(groupMeans(SSVEPTwenyByPP, 2, FrCorrByPP, 'dim'), 4));


colors_1 = [.85, .325, .098;   % Color 1
            .466, .674, .188]; % Color 2

figure();
colors = colors_1; set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');

h = errorBarPlot(sq(permute (SSVEPTwenyByPPByHZCorr, [3 1 2])), 'area', 1, 'alpha', 0.1, 'xaxisvalues', respWindows);


xline(0)
x = [-1200, -900, -900, -1200]; % Define the x-coordinates
y = get(gca, 'YLim'); % Get current y-axis limits
patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% Add a legend to the plot
% Labels correspond to different error and correct conditions, with their locations specified
legend('', 'targeted fr', '', 'non targeted fr', 'Location', 'best');

% Set the title of the plot to describe its content
title('for 20 HZ');


%% resp

SSVEPTwenyByPPByHZCorr = squeeze(nanmean(groupMeans(SSVEPTwenyREspByPP, 2, RFrCorrByPP, 'dim'), 4));


colors_1 = [.85, .325, .098;   % Color 1
            .466, .674, .188]; % Color 2

figure();
colors = colors_1; set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');

h = errorBarPlot(sq(permute (SSVEPTwenyByPPByHZCorr, [3 1 2])), 'area', 1, 'alpha', 0.1, 'xaxisvalues', respWindows);


xline(0)

% Add a legend to the plot
% Labels correspond to different error and correct conditions, with their locations specified
legend('', 'targeted fr', '', 'non targeted fr', 'Location', 'best');

% Set the title of the plot to describe its content
title('for 20 HZ');
%% ploting plots

c=1;
fr_and_e_gen = NaN(size (freq_correct_gen));

for i = 1:2
    for j = 0:1
            a = find ( freq_correct_gen == i & e_or_c_gen == j);
            fr_and_e_gen(a) = c;
            c = c+1;
  
    end
end

FrCorrByPPByCorr = processAndExpand(fr_and_e_gen, pp1, dimSize);

SSVEPTwenyByPPByHZCorr = squeeze(nanmean(groupMeans(SSVEPTwenyByPP, 2, FrCorrByPPByCorr, 'dim'), 4));


colors_1 = [0.4974,0.7088,0.9559; 0.4974,0.7088,0.9559;
     0.5098, 0.3608, 0.7882;0.5098, 0.3608, 0.7882;]; % Color 2

figure();
colors = colors_1; set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');


h = errorBarPlot(sq(permute (SSVEPTwenyByPPByHZCorr, [3 1 2])), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
h{1}.LineStyle = "--"; h{3}.LineStyle = "--";

xline(0)
x = [-1200, -900, -900, -1200]; % Define the x-coordinates
y = get(gca, 'YLim'); % Get current y-axis limits
patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% Add a legend to the plot
% Labels correspond to different error and correct conditions, with their locations specified
legend('', 'targeted fr error', '', ' targeted fr correct', '', 'non-targeted fr error', '', 'non-targeted fr correct', 'Location', 'best');

% Set the title of the plot to describe its content
title('for 20 HZ');

%% ploting plots

c=1;
fr_and_e_ch_gen = NaN(size (freq_ch_gen));

for i = 1:2
    for j = 0:1
            a = find ( freq_ch_gen == i & e_or_c_gen == j);
            fr_and_e_ch_gen(a) = c;
            c = c+1;
  
    end
end

FrCorrByPPByCorr = processAndExpand(fr_and_e_ch_gen, pp1, dimSize);

SSVEPTwenyByPPByHZCorr = squeeze(nanmean(groupMeans(SSVEPTwenyByPP, 2, FrCorrByPPByCorr, 'dim'), 4));


colors_1 = [0.4974,0.7088,0.9559; 0.4974,0.7088,0.9559;
     0.5098, 0.3608, 0.7882;0.5098, 0.3608, 0.7882;]; % Color 2

figure();
colors = colors_1; set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');


h = errorBarPlot(sq(permute (SSVEPTwenyByPPByHZCorr, [3 1 2])), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
h{1}.LineStyle = "--"; h{3}.LineStyle = "--";

xline(0)
x = [-1200, -900, -900, -1200]; % Define the x-coordinates
y = get(gca, 'YLim'); % Get current y-axis limits
patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% Add a legend to the plot
% Labels correspond to different error and correct conditions, with their locations specified
legend('', 'chosen fr error', '', ' chosen fr correct', '', 'non-chosen fr error', '', 'non-chosen fr correct', 'Location', 'best');

% Set the title of the plot to describe its content
title('for 20 HZ');

%% ploting plots

c=1;
fr_and_e_gen = NaN(size (freq_correct_gen));

for i = 1:2
    for j = 0:1
            a = find ( freq_correct_gen == i & e_or_c_gen == j);
            fr_and_e_gen(a) = c;
            c = c+1;
  
    end
end

FrCorrByPPByCorr = processAndExpand(fr_and_e_gen, pp1, dimSize);




SSVEPTwenyByPPByHZCorr = squeeze(nanmean(groupMeans(SSVEPTwenyREspByPP, 2, FrCorrByPPByCorr, 'dim'), 4));


colors_1 = [0.4974,0.7088,0.9559; 0.4974,0.7088,0.9559;
     0.5098, 0.3608, 0.7882;0.5098, 0.3608, 0.7882;]; % Color 2

figure();
colors = colors_1; set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');


h = errorBarPlot(sq(permute (SSVEPTwenyByPPByHZCorr, [3 1 2])), 'area', 1, 'alpha', 0.1, 'xaxisvalues', respWindows);
h{1}.LineStyle = "--"; h{3}.LineStyle = "--";

xline(0)
x = [-1200, -900, -900, -1200]; % Define the x-coordinates
y = get(gca, 'YLim'); % Get current y-axis limits
patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% Add a legend to the plot
% Labels correspond to different error and correct conditions, with their locations specified
legend('', 'target fr error', '', ' target fr correct', '', 'non-target fr error', '', 'non-target fr correct', 'Location', 'best');

% Set the title of the plot to describe its content
title('for 20 HZ');

%% 

SSVEPTwenyFByPPByHZCorr = squeeze(nanmean(groupMeans(SSVEPTwenyFiveByPP, 2, FrCorrByPP, 'dim'), 4));


colors_1 = [.85, .325, .098;   % Color 1
            .466, .674, .188]; % Color 2

figure();
colors = colors_1; set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');

h = errorBarPlot(sq(permute (SSVEPTwenyFByPPByHZCorr, [3 1 2])), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);


xline(0)
x = [-1200, -900, -900, -1200]; % Define the x-coordinates
y = get(gca, 'YLim'); % Get current y-axis limits
patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% Add a legend to the plot
% Labels correspond to different error and correct conditions, with their locations specified
legend('', ' non targeted fr', '', ' targeted fr', 'Location', 'best');

% Set the title of the plot to describe its content
title('for 25 HZ');

%% 
SSVEPTwenyByPPByHZCorr = squeeze(nanmean(groupMeans(SSVEPTwenyFiveByPP, 2, FrCorrByPPByCorr, 'dim'), 4));


colors_1 = [0.4974,0.7088,0.9559; 0.4974,0.7088,0.9559;
     0.5098, 0.3608, 0.7882;0.5098, 0.3608, 0.7882;]; % Color 2

figure();
colors = colors_1; set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');

h = errorBarPlot(sq(permute (SSVEPTwenyByPPByHZCorr, [3 1 2])), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
h{1}.LineStyle = "--"; h{3}.LineStyle = "--";

xline(0)
x = [-1200, -900, -900, -1200]; % Define the x-coordinates
y = get(gca, 'YLim'); % Get current y-axis limits
patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% Add a legend to the plot
% Labels correspond to different error and correct conditions, with their locations specified
legend('', 'non -targeted fr error', '', ' non-targeted fr correct', '', 'targeted fr error', '', 'targeted fr correct', 'Location', 'best');

% Set the title of the plot to describe its content
title('for 25 HZ');
%%

SSVEPTargetByPPByHZCorr = squeeze(nanmean(groupMeans(SSVEP_targetByPP, 2, FrCorrByPP, 'dim'), 4));


colors_1 = [.85, .325, .098;   % Color 1
            .466, .674, .188]; % Color 2

figure();
colors = colors_1; set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');

h = errorBarPlot(sq(permute (SSVEPTargetByPPByHZCorr, [3 1 2])), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);


xline(0)
x = [-1200, -900, -900, -1200]; % Define the x-coordinates
y = get(gca, 'YLim'); % Get current y-axis limits
patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% Add a legend to the plot
% Labels correspond to different error and correct conditions, with their locations specified
legend('', '20', '', '25', 'Location', 'best');

% Set the title of the plot to describe its content
title('target');
%%
SSVEPnonTargetByPPByHZCorr = squeeze(nanmean(groupMeans(SSVEP_nontargetByPP, 2, FrCorrByPP, 'dim'), 4));


colors_1 = [.85, .325, .098;   % Color 1
            .466, .674, .188]; % Color 2

figure();
colors = colors_1; set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');

h = errorBarPlot(sq(permute (SSVEPnonTargetByPPByHZCorr, [3 1 2])), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);


xline(0)
x = [-1200, -900, -900, -1200]; % Define the x-coordinates
y = get(gca, 'YLim'); % Get current y-axis limits
patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% Add a legend to the plot
% Labels correspond to different error and correct conditions, with their locations specified
legend('', '25', '', '20', 'Location', 'best');

% Set the title of the plot to describe its content
title('nontarget');

%% cue
TargatedByPPByCue = squeeze(nanmean(groupMeans(SSVEP_targetByPP, 2, StCueByPP, 'dim'), 4));
NonTargatedByPPByCue = squeeze(nanmean(groupMeans(SSVEP_nontargetByPP, 2, StCueByPP, 'dim'), 4));
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
h = errorBarPlot(movmean(sq(permute (TargatedByPPByCue, [3 1 2])), 5,2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
hold on
e = errorBarPlot(movmean(sq(permute (NonTargatedByPPByCue, [3 1 2])), 5,2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
hold off

e{1}.LineStyle = "--"; e{2}.LineStyle = "--"; e{3}.LineStyle = "--";
xline(0)
yline(0)

% Add a legend to the plot
% Labels correspond to different error and correct conditions, with their locations specified
legend('', 'invalid ', '', 'neutral ','', 'valid ',   'Location', 'best');

% Set the title of the plot to describe its content
title('cue (resp locked)');



%% 
DiffByPPByCue = squeeze(nanmean(groupMeans(SSVEPdiffByPP, 2, StCueByPP, 'dim'), 4));

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
h = errorBarPlot(movmean(sq(permute (DiffByPPByCue, [3 1 2])), 5,2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);

xline(0)
yline(0)
x = [-1200, -900, -900, -1200]; % Define the x-coordinates
y = get(gca, 'YLim'); % Get current y-axis limits
patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% Add a legend to the plot
% Labels correspond to different error and correct conditions, with their locations specified
legend('', 'invalid ', '', 'neutral ','', 'valid ',   'Location', 'best');

% Set the title of the plot to describe its content
title('cue for diff');


%% error/correct
TargetedByPPByAcc = squeeze(nanmean(groupMeans(SSVEP_targetByPP, 2, StAccByPP, 'dim'), 4));
NonTargetedByPPByAcc = squeeze(nanmean(groupMeans(SSVEP_nontargetByPP, 2, StAccByPP, 'dim'), 4));

colors_1 = [.85, .325, .098;   % Color 1
            .466, .674, .188]; % Color 2

figure();

% Define the colors for the plot using colors_5 matrix
colors = colors_1; 

% Set the color order for the plot and specify that new plots replace existing ones in the same axes
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');

% Plot the error bars using the data cppByPPByCueE

% The data is permuted into [trials, participants, cues] for plotting
% 'area' specifies the plotting style, 'alpha' sets the transparency, and 'xaxisvalues' provides x-axis labels
h = errorBarPlot(movmean(sq(permute (TargetedByPPByAcc, [3 1 2])), 5, 2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
hold on
e = errorBarPlot(movmean(sq(permute (NonTargetedByPPByAcc, [3 1 2])), 5,2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
hold off

e{1}.LineStyle = "--"; e{2}.LineStyle = "--"; 
xline(0)
yline(0)

x = [-1200, -900, -900, -1200]; % Define the x-coordinates
y = get(gca, 'YLim'); % Get current y-axis limits
patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

legend('', 'error', '', 'correct', 'Location', 'best');

% Set the title of the plot to describe its content
title('error/correct (stim locked)');

%% certainty
TargetedByPPByCert = squeeze(nanmean(groupMeans(SSVEP_targetByPP, 2, StCertByPP, 'dim'), 4));
NonTargetedByPPByCert = squeeze(nanmean(groupMeans(SSVEP_nontargetByPP, 2, StCertByPP, 'dim'), 4));
colors_1 = [ 0.4974,0.7088,0.9559;
     0.5098, 0.3608, 0.7882;
     0.7804, 0.1412, 0.4000];

figure();

% Define the colors for the plot using colors_5 matrix
colors = colors_1; 

% Set the color order for the plot and specify that new plots replace existing ones in the same axes
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');

h = errorBarPlot(movmean(sq(permute (TargetedByPPByCert(:,:,:), [3 1 2])), 5, 2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
hold on
e = errorBarPlot(movmean(sq(permute (NonTargetedByPPByCert, [3 1 2])), 5,2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
hold off

e{1}.LineStyle = "--"; e{2}.LineStyle = "--"; e{3}.LineStyle = "--"; 
xline(0)
yline(0)

x = [-1200, -900, -900, -1200]; % Define the x-coordinates
y = get(gca, 'YLim'); % Get current y-axis limits
patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

legend('', 'maybe ', '', 'probably  ', '', 'certain ',...
   'Location', 'best');

% Set the title of the plot to describe its content
title('certainty (stim locked)');


%% certainty
DiffByPPByCert = squeeze(nanmean(groupMeans(SSVEPdiffByPP, 2, StCertByPP, 'dim'), 4));

colors_1 = [ 0.4974,0.7088,0.9559;
     0.5098, 0.3608, 0.7882;
     0.7804, 0.1412, 0.4000];

figure();

% Define the colors for the plot using colors_5 matrix
colors = colors_1; 

% Set the color order for the plot and specify that new plots replace existing ones in the same axes
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');

h = errorBarPlot(movmean(sq(permute (DiffByPPByCert(:,:,:), [3 1 2])), 5, 2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);

xline(0)
yline(0)

x = [-1200, -900, -900, -1200]; % Define the x-coordinates
y = get(gca, 'YLim'); % Get current y-axis limits
patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

legend('', 'maybe ', '', 'probably  ', '', 'certain ',...
   'Location', 'best');

% Set the title of the plot to describe its content
title('certainty (stim locked)');
%% block
TargetedByPPBySesion = squeeze(nanmean(groupMeans(SSVEP_targetByPP, 2, SessionByPP, 'dim'), 4));
NonTargetedByPPBySesion = squeeze(nanmean(groupMeans(SSVEP_nontargetByPP, 2, SessionByPP, 'dim'), 4));
colors_3 = [0.925, 0.6625, 0.549;  % Color 1
            0.85, 0.325, 0.098;    % Color 2
            0.425, 0.1625, 0.049
             0.5, 0.724, 0.8705;    % Color 4
            0, 0.447, 0.741;       % Color 5
            0, 0.2235, 0.3705; ]; % Color 9
figure();

% Define the colors for the plot using colors_5 matrix
colors = colors_3; 

% Set the color order for the plot and specify that new plots replace existing ones in the same axes
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');

h = errorBarPlot(movmean(sq(permute (TargetedByPPBySesion(:,:,:), [3 1 2])), 5, 2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
hold on
e = errorBarPlot(movmean(sq(permute (NonTargetedByPPBySesion, [3 1 2])), 5,2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
hold off

e{1}.LineStyle = "--"; e{2}.LineStyle = "--"; e{3}.LineStyle = "--"; 
xline(0)
yline(0)

x = [-1200, -900, -900, -1200]; % Define the x-coordinates
y = get(gca, 'YLim'); % Get current y-axis limits
patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

legend('', '1 ', '', '2  ', '', '3 ',...
   'Location', 'best');

% Set the title of the plot to describe its content
title('certainty (stim locked)');

%% 1 block+cue
cue_m_gen1 = cue_m_gen;
cue_m_gen1(session_gen ~= 1) = NaN;

StCueByPPBySess1 = processAndExpand(cue_m_gen1, pp1, dimSize);

TargatedByPPByCueSess = squeeze(nanmean(groupMeans(SSVEP_targetByPP, 2, StCueByPPBySess1, 'dim'), 4));
NonTargatedByPPByCueSess = squeeze(nanmean(groupMeans(SSVEP_nontargetByPP, 2, StCueByPPBySess1, 'dim'), 4));
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
h = errorBarPlot(movmean(sq(permute (TargatedByPPByCueSess, [3 1 2])), 5,2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
hold on
e = errorBarPlot(movmean(sq(permute (NonTargatedByPPByCueSess, [3 1 2])), 5,2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
hold off

e{1}.LineStyle = "--"; e{2}.LineStyle = "--"; e{3}.LineStyle = "--";
xline(0)
yline(0)
ylim ([-0.1 0.7])
% Add a legend to the plot
% Labels correspond to different error and correct conditions, with their locations specified
legend('', 'invalid ', '', 'neutral ','', 'valid ',   'Location', 'best');

% Set the title of the plot to describe its content
title('cue (1 session)');

%% 2 block+cue
cue_m_gen1 = cue_m_gen;
cue_m_gen1(session_gen ~= 2) = NaN;

StCueByPPBySess1 = processAndExpand(cue_m_gen1, pp1, dimSize);

TargatedByPPByCueSess = squeeze(nanmean(groupMeans(SSVEP_targetByPP, 2, StCueByPPBySess1, 'dim'), 4));
NonTargatedByPPByCueSess = squeeze(nanmean(groupMeans(SSVEP_nontargetByPP, 2, StCueByPPBySess1, 'dim'), 4));
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
h = errorBarPlot(movmean(sq(permute (TargatedByPPByCueSess, [3 1 2])), 5,2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
hold on
e = errorBarPlot(movmean(sq(permute (NonTargatedByPPByCueSess, [3 1 2])), 5,2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
hold off

e{1}.LineStyle = "--"; e{2}.LineStyle = "--"; e{3}.LineStyle = "--";
xline(0)
yline(0)
ylim ([-0.1 0.7])
% Add a legend to the plot
% Labels correspond to different error and correct conditions, with their locations specified
legend('', 'invalid ', '', 'neutral ','', 'valid ',   'Location', 'best');

% Set the title of the plot to describe its content
title('cue (2 session)');


%% 3 block+cue
cue_m_gen1 = cue_m_gen;
cue_m_gen1(session_gen ~= 3) = NaN;

StCueByPPBySess1 = processAndExpand(cue_m_gen1, pp1, dimSize);

TargatedByPPByCueSess = squeeze(nanmean(groupMeans(SSVEP_targetByPP, 2, StCueByPPBySess1, 'dim'), 4));
NonTargatedByPPByCueSess = squeeze(nanmean(groupMeans(SSVEP_nontargetByPP, 2, StCueByPPBySess1, 'dim'), 4));
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
h = errorBarPlot(movmean(sq(permute (TargatedByPPByCueSess, [3 1 2])), 5,2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
hold on
e = errorBarPlot(movmean(sq(permute (NonTargatedByPPByCueSess, [3 1 2])), 5,2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);
hold off

e{1}.LineStyle = "--"; e{2}.LineStyle = "--"; e{3}.LineStyle = "--";
xline(0)
yline(0)
ylim ([-0.1 0.7])
% Add a legend to the plot
% Labels correspond to different error and correct conditions, with their locations specified
legend('', 'invalid ', '', 'neutral ','', 'valid ',   'Location', 'best');

% Set the title of the plot to describe its content
title('cue (3 session)');


%% certainty 1 block

cue_m_gen1 = cert_gen;
cue_m_gen1(session_gen ~= 1) = NaN;


StCertByPPBySess1 = processAndExpand(cue_m_gen1, pp1, dimSize);


DiffByPPByCue = squeeze(nanmean(groupMeans(SSVEPdiffByPP, 2, StCertByPPBySess1, 'dim'), 4));

colors_1 = [ 0.4974,0.7088,0.9559;
     0.5098, 0.3608, 0.7882;
     0.7804, 0.1412, 0.4000];

figure();

% Define the colors for the plot using colors_5 matrix
colors = colors_1; 

% Set the color order for the plot and specify that new plots replace existing ones in the same axes
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');

h = errorBarPlot(movmean(sq(permute (DiffByPPByCue(:,:,:), [3 1 2])), 5, 2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);

xline(0)
yline(0)

ylim ([-0.1 0.7])


x = [-1200, -900, -900, -1200]; % Define the x-coordinates
y = get(gca, 'YLim'); % Get current y-axis limits
patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

legend('', 'maybe ', '', 'probably  ', '', 'certain ',...
   'Location', 'best');

% Set the title of the plot to describe its content
title('certainty (1 session)');


%% 


cue_m_gen1 = cue_m_gen;
cue_m_gen1(session_gen ~= 1) = NaN;


StCertByPPBySess1 = processAndExpand(cue_m_gen1, pp1, dimSize);

DiffByPPByCue = squeeze(nanmean(groupMeans(SSVEPdiffByPP, 2, StCertByPPBySess1, 'dim'), 4));

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
h = errorBarPlot(movmean(sq(permute (DiffByPPByCue, [3 1 2])), 5,2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);

xline(0)
yline(0)
x = [-1200, -900, -900, -1200]; % Define the x-coordinates
y = get(gca, 'YLim'); % Get current y-axis limits
patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% Add a legend to the plot
% Labels correspond to different error and correct conditions, with their locations specified
legend('', 'invalid ', '', 'neutral ','', 'valid ',   'Location', 'best');

% Set the title of the plot to describe its content
title('cue for diff');
%% certainty 2 block

cue_m_gen1 = cert_gen;
cue_m_gen1(session_gen ~= 2) = NaN;


StCertByPPBySess1 = processAndExpand(cue_m_gen1, pp1, dimSize);


DiffByPPByCue = squeeze(nanmean(groupMeans(SSVEPdiffByPP, 2, StCertByPPBySess1, 'dim'), 4));
colors_1 = [ 0.4974,0.7088,0.9559;
     0.5098, 0.3608, 0.7882;
     0.7804, 0.1412, 0.4000];

figure();

% Define the colors for the plot using colors_5 matrix
colors = colors_1; 

% Set the color order for the plot and specify that new plots replace existing ones in the same axes
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');

h = errorBarPlot(movmean(sq(permute (DiffByPPByCue(:,:,:), [3 1 2])), 5, 2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);

xline(0)
yline(0)


x = [-1200, -900, -900, -1200]; % Define the x-coordinates
y = get(gca, 'YLim'); % Get current y-axis limits
patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

legend('', 'maybe ', '', 'probably  ', '', 'certain ',...
   'Location', 'best');

% Set the title of the plot to describe its content
title('certainty (2 session)');


%% 


cue_m_gen1 = cue_m_gen;
cue_m_gen1(session_gen ~= 2) = NaN;


StCertByPPBySess1 = processAndExpand(cue_m_gen1, pp1, dimSize);

DiffByPPByCue = squeeze(nanmean(groupMeans(SSVEPdiffByPP, 2, StCertByPPBySess1, 'dim'), 4));

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
h = errorBarPlot(movmean(sq(permute (DiffByPPByCue, [3 1 2])), 5,2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);

xline(0)
yline(0)
x = [-1200, -900, -900, -1200]; % Define the x-coordinates
y = get(gca, 'YLim'); % Get current y-axis limits
patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% Add a legend to the plot
% Labels correspond to different error and correct conditions, with their locations specified
legend('', 'invalid ', '', 'neutral ','', 'valid ',   'Location', 'best');

% Set the title of the plot to describe its content
title('cue for diff');
%% certainty 3 block

cue_m_gen1 = cert_gen;
cue_m_gen1(session_gen ~= 3) = NaN;


StCertByPPBySess1 = processAndExpand(cue_m_gen1, pp1, dimSize);


DiffByPPByCue = squeeze(nanmean(groupMeans(SSVEPdiffByPP, 2, StCertByPPBySess1, 'dim'), 4));

colors_1 = [ 0.4974,0.7088,0.9559;
     0.5098, 0.3608, 0.7882;
     0.7804, 0.1412, 0.4000];

figure();

% Define the colors for the plot using colors_5 matrix
colors = colors_1; 

% Set the color order for the plot and specify that new plots replace existing ones in the same axes
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');

h = errorBarPlot(movmean(sq(permute (DiffByPPByCue(:,:,:), [3 1 2])), 5, 2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);

xline(0)
yline(0)


x = [-1200, -900, -900, -1200]; % Define the x-coordinates
y = get(gca, 'YLim'); % Get current y-axis limits
patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

legend('', 'maybe ', '', 'probably  ', '', 'certain ',...
   'Location', 'best');

% Set the title of the plot to describe its content
title('certainty (3 session)');

%% 


cue_m_gen1 = cue_m_gen;
cue_m_gen1(session_gen ~= 3) = NaN;


StCertByPPBySess1 = processAndExpand(cue_m_gen1, pp1, dimSize);

DiffByPPByCue = squeeze(nanmean(groupMeans(SSVEPdiffByPP, 2, StCertByPPBySess1, 'dim'), 4));

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
h = errorBarPlot(movmean(sq(permute (DiffByPPByCue, [3 1 2])), 5,2), 'area', 1, 'alpha', 0.1, 'xaxisvalues', stimWindows);

xline(0)
yline(0)
x = [-1200, -900, -900, -1200]; % Define the x-coordinates
y = get(gca, 'YLim'); % Get current y-axis limits
patch(x, [y(1) y(1) y(2) y(2)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% Add a legend to the plot
% Labels correspond to different error and correct conditions, with their locations specified
legend('', 'invalid ', '', 'neutral ','', 'valid ',   'Location', 'best');

% Set the title of the plot to describe its content
title('cue for diff');
