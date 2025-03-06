function my_CheckForMismatch(folders)


clc


if ~exist('pop_biosig','file')
    eeglab; close all;
end

close all;

help Edf2Mat

%% loop for getting timing and responses from matlab, eeg and eye-tracker
%  for i = 1:numel(pid)
% folderNames = [pid{i}];
folder1 = fullfile('D:\cue_task\analysis\Data\', folders{1});
% 
%  end
%folder1 = ['C:\Users\monakhvd\OneDrive - Trinity College Dublin\Data\P09\1']

mat_files = dir(fullfile(folder1, '*2024*.mat'));
bdf_files = dir(fullfile(folder1, '*.bdf'));
edf_files = dir(fullfile(folder1, '*.edf'));


%matrix for timing and trigs from matlab, eeg and eye-tracker

time_and_trigs = struct(...
    'eeg_t', [], 'mat_t', [], ...
    'eye_t', [], 'mat_trigs', [], ...
    'eeg_trigs', [], 'eye_trigs', [] ...
);
sorted_table = [];
%% loop through all blocks

for i = 1 :6
%loading matlab
    load (fullfile(folder1, mat_files(i).name));
    PTBT = [PTBeventT]-PTBeventT(1); %storring timing with correction for the 1 trigger

%loading eeg and getting trig/timing
    time_and_trigs.mat_t = [time_and_trigs.mat_t, PTBT]; % getting timing from matlab file
    time_and_trigs.mat_trigs = [time_and_trigs.mat_trigs, PTBevent]; %getting triggers 
    f = time_and_trigs.mat_trigs (2:8)'; f2 = time_and_trigs.mat_t(2:8)'; f2 = f2 - f2(1); f2= round (f2*1000);
checkiing_time_for_t = [f, f2]; 

     EEG = pop_biosig (fullfile(folder1, bdf_files(i).name));
    disp(['Amount of triggers in EEG for block ', num2str(i), ': ', num2str(length ([EEG.event.urevent]))]);
    disp(['Amount of triggers in MATLAB for block ', num2str(i), ': ', num2str(length (PTBevent))]); 
      %getting triggers from EEG
nEvents = length(PTBevent)-1;
trigs = [EEG.event.edftype]; % get number triggers
a = sum(trigs==1);

% as some triggers arre not saved by in trigs, we need to insert them

        if length(trigs) < nEvents % if fewer triggers found than expected
     
            emptyInds = find(cellfun(@isempty, {EEG.event.edftype})); % find missing ones

            for j = 1:length(emptyInds) % convert them from the string in .type
                EEG.event(emptyInds(j)).edftype = str2double(EEG.event(emptyInds(j)).type);
            end
        
            trigs = [EEG.event.edftype]; % get them again

            if length(trigs) < nEvents % if there are still missing ones

                trigs = {EEG.event.type}; % this is sometimes 'condition X' or just 'X' (where X is a string number)
                if iscellstr(trigs)
                    % get triggers as numbers only
                    trigsNum = NaN(1,nEvents);
                    for k = 1:nEvents
                        trigsNum(k) = str2double(trigs{k}(regexp(trigs{k}, '\d')));
                    end

                    trigs = trigsNum; % replace
                elseif all(isnumeric([trigs{:}]))
                    trigs = [trigs{:}]; %
                else

                    keyboard;
                    % error('trigs not converted to numbers');
                end
            end
        end

    [c,k] = find (trigs > 40); % finding odd triggers and where they are
    disp([c; k]);
    disp (trigs (k))
    % keyboard


lat = [EEG.event.latency] ./ 1024; % sampling down eeg timing
time_and_trigs.eeg_trigs = [time_and_trigs.eeg_trigs, trigs]; % saving triggers from eeg +getting rid of odd triggers
lat = lat - lat(1); % storring timing with correction for the 1 trigger
time_and_trigs.eeg_t = [time_and_trigs.eeg_t, lat]; %storring all timing in one array



   edf1 = Edf2Mat(fullfile(folder1, edf_files(i).name));
   a = edf1.Events.Messages.time; %getting time out of eye-tracker
   b = edf1.Events.Messages.info; %getting trigs out of eye-tracker
   trigs_eye = [];
   number = [];

   for j = 1 : (length (b))
    c = cell2mat(b(j));
    match = regexp(c, 'TRIG\d+', 'match', 'once');

    if ~isempty(match)
    firstNumber = str2double(regexp(match, '\d+', 'match', 'once'));

    % if firstNumber >= 1 && firstNumber <= 40 
        trigs_eye = [trigs_eye, firstNumber];
        number = [number, j];
    %end

    end
end

disp(['Amount of triggers in EYE-TRAKER for block ', num2str(i), ': ', num2str(length (trigs_eye))]);

%getting timing from eye-tracker
t_eye = [];
for r = 1 : length (number)
    timimg = a (number(r));
    t_eye = [t_eye, timimg];
end
    t_eye = t_eye - edf1.Events.Start.time(1);
    t_eye = t_eye/1000;

    time_and_trigs.eye_t = [time_and_trigs.eye_t, t_eye];
    time_and_trigs.eye_trigs = [time_and_trigs.eye_trigs, trigs_eye];
   
    %checking if there was mismatch
   % Ensure unmatched_variable is defined
unmatched_variable = find(trigs(1:length(trigs_eye)) ~= trigs_eye);

if ~isempty(unmatched_variable)  % Check if unmatched_variable is non-empty
    for z = 1:length(unmatched_variable)
    %for z = 1:1000
        % Fix the condition to check for unmatched triggers
        if (trigs(unmatched_variable(z)) == trigs_eye(unmatched_variable(z))) || ...
           (z > 1 && trigs(unmatched_variable(z)) == trigs_eye(unmatched_variable(z-1)))
            continue;  % Skip the iteration if condition is met
        else
            l = trigs(unmatched_variable(z)-6:unmatched_variable(z)+3);
            k = trigs_eye(unmatched_variable(z)-6:unmatched_variable(z)+3);
            disp(['EEG trig: ', num2str(l)]);
            disp(['ET trig: ', num2str(k)]);
            fprintf('\nnumber: %.2f\n', unmatched_variable(z));
            
            % Extract participant name
            p_name = extractBefore(mat_files(i).name, '_Test');
            
           % Determine the write action
          if trigs(unmatched_variable(z)) - trigs_eye(unmatched_variable(z)) == 32 || ...
            trigs_eye(unmatched_variable(z)) - trigs(unmatched_variable(z) - 1) == 32
    writie = 'change trig';
          elseif length(trigs) > length(trigs_eye) && ...
        trigs(unmatched_variable(z)) ~= trigs_eye(unmatched_variable(z)) && ...
        trigs(unmatched_variable(z)) ~= trigs_eye(unmatched_variable(z) - 1)
    writie = 'delete trig';
         
end
            % 
            % Update sorted_table
            sorted_table = [sorted_table; {p_name}, num2cell(unmatched_variable(z)) ...
                num2cell(trigs(unmatched_variable(z))), ...
                num2cell(trigs_eye(unmatched_variable(z))), {writie}]; %
        end
    end
end
end

if ~isempty(sorted_table)
what_should_be_del = []; %Calculate the percentage of responses relative to the total number of occurrences
what_should_be_del = array2table ([sorted_table], 'VariableNames',{ 'block', 'ind',  'trig', ...
    'eye trig', 'to do'}) %table with names for them
%keyboard
end


% 
% if ~isempty(sorted_table)
% save (fullfile(folder1, [p_name '_trig_mismatch']), 'what_should_be_del')
% end