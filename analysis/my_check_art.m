clear all
clc


fold = 'C:\Users\monakhvd\OneDrive - Trinity College Dublin\Data\Combine';

f = what(fold);
wh = 'P17';
% files = f.mat(cellRegexpi(f.mat, '15.23')>0);
pattern = ['.*', wh, '.*flagged.*'];
files = f.mat(cellRegexpi(f.mat, pattern)>0);
erp_gen = [];isGood_gen = []; RT = [];

for i = 1:length (files)
    
    load (fullfile(fold, cell2mat (files (i))));

   
    isGood = isGood(1:864);
    RT = [RT; RS];

    %cleaning ERP
isGood_gen = [isGood_gen;isGood ];

end 
length (find (isGood_gen == 1))
% folder = 'C:\Users\monakhvd\OneDrive - Trinity College Dublin\Data\General';
% fprintf('savig...');
% save (fullfile(folder, [wh '_whole']),'-v7.3',...
%         'isGood_gen', 'erp_gen', 'RT');
% fprintf('end');

