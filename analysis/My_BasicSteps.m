%%%% general script of my plots

% 
%% checking 1 session at a time
clear all
folder1 = 'D:\cue_task\analysis\Data\P21\beh'; % changelast bit of path
wh = 'P21_1'; % change to wanted session number
my_BehCheck

%% making one file for cpp and flagged
my_combining %for combining, need to check which participant is there

%% checking 3 sessions together (only beh)
close all; clear all; 
% 
% %choose participant 
 wh = ['P23'];

 my_BehThreeSessions

%% checking 3 sessions together of erp and saving cpp
myERPfull

%% mean over all participants
clear all
sub = {'P01' 'P02' 'P03' 'P04'  'P05' 'P06' 'P07' 'P08' ...
    'P10'  'P13'  'P14' 'P15' 'P16' 'P17' 'P18' 'P19' 'P20' 'P21' 'P22'}; %delete if needed check if it's commented in actual code

my_BehMeanOverAll

%% conditional plot for accuracy over time
% sub = {'P01' 'P02' 'P03' 'P04'  'P05' 'P06' 'P07' 'P08' ...
%     'P10'  'P13'  'P14' 'P15' 'P16' 'P17' 'P18' 'P19' 'P20' 'P21'};%'P20'check if it's commented in actual code

MyPlotingAccuacyAsRTfuncs

%% building general plots of CPP

MyCPP % for resp
MyCPP_stim %for stim

% sub = {'P01' 'P02' 'P03' 'P04'  'P05' 'P06' 'P07' 'P08' ...
%     'P10'  'P13'  'P14' 'P15' 'P16' 'P17' 'P18' 'P19' 'P20' 'P21'};%'P20'check if it's commented in actual code
MyLoopThroughEachCPP % to plot each person separatly
%% conditional plot for CPP over time
% sub = {'P01' 'P02' 'P03' 'P04'  'P05' 'P06' 'P07' 'P08' ...
%     'P10'  'P13'  'P14' 'P15' 'P16' 'P17' 'P18' 'P19' 'P20' 'P21'};%'P20' check if it's commented in actual code
MyPlotingCPPasRTfuncs


%% mu beta

my_beta_plotting

%% SSVEP

%topoplot for ssvep
my_twenty_twentyFive

my_SSVEP_save

my_SSVEP_plot

%% RMANOVA and modeling

my_modeling