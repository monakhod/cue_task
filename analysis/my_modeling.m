clc
clear all

folder = 'G:\cue_task\analysis\Data\General';
folder1 = 'G:\cue_task\analysis\Data\means';


%% loading each participant

behTab_gen = [];
mat_files = dir(fullfile(folder, '*sorted_t.mat'));

for i = 1:length (mat_files)
    fprintf('Loading files: %s \n ', mat_files(i).name);
    load (fullfile(folder,  (mat_files(i).name)));
    behTab_gen = [behTab_gen; behTab];
    
end 

load (fullfile(folder1,  'mean_ccp.mat'));
behTab_gen = addvars(behTab_gen, eeg_mean, 'NewVariableNames', 'cpp');


%% repeated measures anovas
% factors
factorNames = {'cue'};
% behTab loaded up
regTab = behTab_gen; %clear behTab;
 % remove all NaN rows
rowsWithNaN = any(isnan(table2array(regTab(:,2))), 2);
regTab(rowsWithNaN, :) = [];
regNames = regTab.Properties.VariableNames;

for i = 1:length(regNames)
    if  i == 2 | i == 6 | i == 8 | i == 9
        continue; % Skip the current iteration for indices 2 and 5
    end
    regTab.(regNames{i}) = nanzscore(regTab.(regNames{i}));
end

logIVs = {'initAcc'} %,'confAcc','CoM'};
ivNames = regNames;
logInds = find(ismember(regNames, logIVs));
ivNames(logInds) = strcat(ivNames(logInds), 'Logistic');

isLogistic = ismember(regNames, ivNames(logInds)); % which are [0 1] for logistic regressions

%% regress by cue
% vars
varNames = {'initAccLogistic','RT','certainty', 'cpp'};
 nVars = length(varNames); % just these
lmeArgs = { {}, {'link','logit','distribution','binomial'} }; % for normal/logistic regression

formulas = '%s ~ 1 + cue + (1 | pp)';


maxInteraction = []; % max order of interactions tested in the RE structure. []=set by the FE terms in formula. or can set to 0 to ignore all RE-interactions
alphaAccept = 0.05; % p-threshold for adding a term

for iV = 1: nVars
    [fits{iV}, formulas1{iV}] = stepUpRECorrel(regTab, sprintf(formulas, varNames{iV}), lmeArgs{ isLogistic(strcmp(regNames, varNames{iV})) + 1}{:}, maxInteraction, alphaAccept);
    % steps up through (correlated) random slope terms, adding ones that improve the model fit, then tests against an uncorrelated version of the best one
    % or you can use stepUpRE() if you prefer to step-wise through the uncorrelated terms and then tests against the correlated one, though I find StepUpRECorrel works better
end

regTables = struct();
fits = cell(nVars,1);
for iV = 1: nVars
    fits{iV} = fitglme(regTab, sprintf(formulas, varNames{iV}), lmeArgs{ isLogistic(strcmp(regNames, varNames{iV})) + 1}{:});
end
regTables.p = array2table(sq(nancat(cellfun(@(x) x.Coefficients.pValue(2:end), fits, 'UniformOutput',0)))','VariableNames',varNames,'RowNames',fits{1}.CoefficientNames(2:end)); 
regTables.beta = array2table(sq(nancat(cellfun(@(x) x.Coefficients.Estimate(2:end), fits, 'UniformOutput',0)))','VariableNames',varNames,'RowNames',fits{1}.CoefficientNames(2:end)); 

%% regress by cert
% vars
varNames = {'initAccLogistic','RT', 'cpp'}; 
 nVars = length(varNames); % just these
lmeArgs = { {}, {'link','logit','distribution','binomial'} }; % for normal/logistic regression

formulas5 = '%s ~ 1 + certainty + (1 | pp) + (-1 + certainty | pp)';

for iV = 1: nVars
    [fits{iV}, formulas1{iV}] = stepUpRECorrel(regTab, sprintf(formulas5, varNames{iV}), lmeArgs{ isLogistic(strcmp(regNames, varNames{iV})) + 1}{:}, maxInteraction, alphaAccept);
    % steps up through (correlated) random slope terms, adding ones that improve the model fit, then tests against an uncorrelated version of the best one
    % or you can use stepUpRE() if you prefer to step-wise through the uncorrelated terms and then tests against the correlated one, though I find StepUpRECorrel works better
end

regTables1 = struct();
fits1 = cell(nVars,1);
for iV = 1: nVars
    fits1{iV} = fitglme(regTab, sprintf(formulas5, varNames{iV}), lmeArgs{ isLogistic(strcmp(regNames, varNames{iV})) + 1}{:});
end
regTables1.p = array2table(sq(nancat(cellfun(@(x) x.Coefficients.pValue(2:end), fits1, 'UniformOutput',0)))','VariableNames',varNames,'RowNames',fits1{1}.CoefficientNames(2:end)); 
regTables1.beta = array2table(sq(nancat(cellfun(@(x) x.Coefficients.Estimate(2:end), fits1, 'UniformOutput',0)))','VariableNames',varNames,'RowNames',fits1{1}.CoefficientNames(2:end)); 


%%  interaction of cue*cert

varNames = {'RT', 'initAccLogistic', 'cpp'}; nVars = length(varNames); % just these
formulas2 = '%s~ 1+ cue*certainty  + (1 | pp)+ (-1 + certainty| pp)';

maxInteraction = []; % max order of interactions tested in the RE structure. []=set by the FE terms in formula. or can set to 0 to ignore all RE-interactions
alphaAccept = 0.05; % p-threshold for adding a term

for iV = 1: nVars
    [fits{iV}, formulas1{iV}] = stepUpRECorrel(regTab, sprintf(formulas2, varNames{iV}), lmeArgs{ isLogistic(strcmp(regNames, varNames{iV})) + 1}{:}, maxInteraction, alphaAccept);
    % steps up through (correlated) random slope terms, adding ones that improve the model fit, then tests against an uncorrelated version of the best one
    % or you can use stepUpRE() if you prefer to step-wise through the uncorrelated terms and then tests against the correlated one, though I find StepUpRECorrel works better
end

regTables2 = struct();

fits2 = cell(nVars,1);
for iV = 1:nVars
    fits2{iV} = fitglme(regTab, sprintf(formulas2, varNames{iV}), lmeArgs{ isLogistic(strcmp(regNames, varNames{iV})) + 1}{:});
end
regTables2.p = array2table(sq(nancat(cellfun(@(x) x.Coefficients.pValue(2:end), fits2, 'UniformOutput',0))),'VariableNames',varNames,'RowNames',fits2{1}.CoefficientNames(2:end)); 
regTables2.beta = array2table(sq(nancat(cellfun(@(x) x.Coefficients.Estimate(2:end), fits2, 'UniformOutput',0))),'VariableNames',varNames,'RowNames',fits2{1}.CoefficientNames(2:end)); 

%%  interaction of cue*initAcc

varNames = {'RT', 'certainty', 'cpp'}; nVars = length(varNames); % just these
formulas2 = '%s~ 1+ cue+initAcc  + (1 | pp)+ (1   +initAcc | pp)';

maxInteraction = []; % max order of interactions tested in the RE structure. []=set by the FE terms in formula. or can set to 0 to ignore all RE-interactions
alphaAccept = 0.05; % p-threshold for adding a term

for iV = 1: nVars
    [fits{iV}, formulas1{iV}] = stepUpRECorrel(regTab, sprintf(formulas2, varNames{iV}), lmeArgs{ isLogistic(strcmp(regNames, varNames{iV})) + 1}{:}, maxInteraction, alphaAccept);
    % steps up through (correlated) random slope terms, adding ones that improve the model fit, then tests against an uncorrelated version of the best one
    % or you can use stepUpRE() if you prefer to step-wise through the uncorrelated terms and then tests against the correlated one, though I find StepUpRECorrel works better
end

regTables3 = struct();

fits3 = cell(nVars,1);
for iV = 1:nVars
    fits3{iV} = fitglme(regTab, sprintf(formulas2, varNames{iV}), lmeArgs{ isLogistic(strcmp(regNames, varNames{iV})) + 1}{:});
end
regTables3.p = array2table(sq(nancat(cellfun(@(x) x.Coefficients.pValue(2:end), fits3, 'UniformOutput',0))),'VariableNames',varNames,'RowNames',fits3{1}.CoefficientNames(2:end)); 
regTables3.beta = array2table(sq(nancat(cellfun(@(x) x.Coefficients.Estimate(2:end), fits3, 'UniformOutput',0))),'VariableNames',varNames,'RowNames',fits3{1}.CoefficientNames(2:end)); 

%%  interaction of cue*RT

varNames = {'initAccLogistic', 'certainty', 'cpp'}; nVars = length(varNames); % just these
formulas2 = '%s~ 1+ cue*RT  + (1 | pp)+ (-1 + RT | pp)';

maxInteraction = []; % max order of interactions tested in the RE structure. []=set by the FE terms in formula. or can set to 0 to ignore all RE-interactions
alphaAccept = 0.05; % p-threshold for adding a term

for iV = 1: nVars
    [fits{iV}, formulas1{iV}] = stepUpRECorrel(regTab, sprintf(formulas2, varNames{iV}), lmeArgs{ isLogistic(strcmp(regNames, varNames{iV})) + 1}{:}, maxInteraction, alphaAccept);
    % steps up through (correlated) random slope terms, adding ones that improve the model fit, then tests against an uncorrelated version of the best one
    % or you can use stepUpRE() if you prefer to step-wise through the uncorrelated terms and then tests against the correlated one, though I find StepUpRECorrel works better
end

regTables4 = struct();

fits4 = cell(nVars,1);
for iV = 1:nVars
    fits4{iV} = fitglme(regTab, sprintf(formulas2, varNames{iV}), lmeArgs{ isLogistic(strcmp(regNames, varNames{iV})) + 1}{:});
end
regTables4.p = array2table(sq(nancat(cellfun(@(x) x.Coefficients.pValue(2:end), fits4, 'UniformOutput',0))),'VariableNames',varNames,'RowNames',fits4{1}.CoefficientNames(2:end)); 
regTables4.beta = array2table(sq(nancat(cellfun(@(x) x.Coefficients.Estimate(2:end), fits4, 'UniformOutput',0))),'VariableNames',varNames,'RowNames',fits4{1}.CoefficientNames(2:end)); 

%%  interaction of cue*RT*initAcc

varNames = {'certainty', 'cpp'}; nVars = length(varNames); % just these
formulas2 = '%s~ 1+ cue*RT*initAcc  + (1 | pp)';

maxInteraction = []; % max order of interactions tested in the RE structure. []=set by the FE terms in formula. or can set to 0 to ignore all RE-interactions
alphaAccept = 0.05; % p-threshold for adding a term

for iV = 1: nVars
    [fits{iV}, formulas1{iV}] = stepUpRECorrel(regTab, sprintf(formulas2, varNames{iV}), lmeArgs{ isLogistic(strcmp(regNames, varNames{iV})) + 1}{:}, maxInteraction, alphaAccept);
    % steps up through (correlated) random slope terms, adding ones that improve the model fit, then tests against an uncorrelated version of the best one
    % or you can use stepUpRE() if you prefer to step-wise through the uncorrelated terms and then tests against the correlated one, though I find StepUpRECorrel works better
end

regTables5 = struct();

fits5 = cell(nVars,1);
for iV = 1:nVars
    fits5{iV} = fitglme(regTab, sprintf(formulas2, varNames{iV}), lmeArgs{ isLogistic(strcmp(regNames, varNames{iV})) + 1}{:});
end
regTables5.p = array2table(sq(nancat(cellfun(@(x) x.Coefficients.pValue(2:end), fits5, 'UniformOutput',0))),'VariableNames',varNames,'RowNames',fits5{1}.CoefficientNames(2:end)); 
regTables5.beta = array2table(sq(nancat(cellfun(@(x) x.Coefficients.Estimate(2:end), fits5, 'UniformOutput',0))),'VariableNames',varNames,'RowNames',fits5{1}.CoefficientNames(2:end)); 

%%  interaction of cue*RT*initAcc

varNames = {'certainty', 'cpp'}; nVars = length(varNames); % just these
formulas2 = '%s~ 1+ initAcc*RT + initAcc*cue + RT*cue + initAcc:RT:cue + (1  +RT  +initAcc | pp)';

maxInteraction = []; % max order of interactions tested in the RE structure. []=set by the FE terms in formula. or can set to 0 to ignore all RE-interactions
alphaAccept = 0.05; % p-threshold for adding a term

for iV = 1: nVars
    [fits{iV}, formulas1{iV}] = stepUpRECorrel(regTab, sprintf(formulas2, varNames{iV}), lmeArgs{ isLogistic(strcmp(regNames, varNames{iV})) + 1}{:}, maxInteraction, alphaAccept);
    % steps up through (correlated) random slope terms, adding ones that improve the model fit, then tests against an uncorrelated version of the best one
    % or you can use stepUpRE() if you prefer to step-wise through the uncorrelated terms and then tests against the correlated one, though I find StepUpRECorrel works better
end

regTables6 = struct();

fits6 = cell(nVars,1);
for iV = 1:nVars
    fits6{iV} = fitglme(regTab, sprintf(formulas2, varNames{iV}), lmeArgs{ isLogistic(strcmp(regNames, varNames{iV})) + 1}{:});
end
regTables6.p = array2table(sq(nancat(cellfun(@(x) x.Coefficients.pValue(2:end), fits6, 'UniformOutput',0))),'VariableNames',varNames,'RowNames',fits6{1}.CoefficientNames(2:end)); 
regTables6.beta = array2table(sq(nancat(cellfun(@(x) x.Coefficients.Estimate(2:end), fits6, 'UniformOutput',0))),'VariableNames',varNames,'RowNames',fits6{1}.CoefficientNames(2:end)); 

%%  interaction of cue*certainty*initAcc

varNames = {'RT', 'cpp'}; nVars = length(varNames); % just these
formulas2 = '%s~ 1+ cue*certainty*initAcc  + (1 | pp)+ (-1 + certainty| pp)';

maxInteraction = []; % max order of interactions tested in the RE structure. []=set by the FE terms in formula. or can set to 0 to ignore all RE-interactions
alphaAccept = 0.05; % p-threshold for adding a term

for iV = 1: nVars
    [fits{iV}, formulas1{iV}] = stepUpRECorrel(regTab, sprintf(formulas2, varNames{iV}), lmeArgs{ isLogistic(strcmp(regNames, varNames{iV})) + 1}{:}, maxInteraction, alphaAccept);
    % steps up through (correlated) random slope terms, adding ones that improve the model fit, then tests against an uncorrelated version of the best one
    % or you can use stepUpRE() if you prefer to step-wise through the uncorrelated terms and then tests against the correlated one, though I find StepUpRECorrel works better
end

regTables7 = struct();

fits7 = cell(nVars,1);
for iV = 1:nVars
    fits7{iV} = fitglme(regTab, sprintf(formulas2, varNames{iV}), lmeArgs{ isLogistic(strcmp(regNames, varNames{iV})) + 1}{:});
end
regTables7.p = array2table(sq(nancat(cellfun(@(x) x.Coefficients.pValue(2:end), fits7, 'UniformOutput',0))),'VariableNames',varNames,'RowNames',fits7{1}.CoefficientNames(2:end)); 
regTables7.beta = array2table(sq(nancat(cellfun(@(x) x.Coefficients.Estimate(2:end), fits7, 'UniformOutput',0))),'VariableNames',varNames,'RowNames',fits7{1}.CoefficientNames(2:end)); 


