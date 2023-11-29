% This one's easy. It sucks in and compiles all the data and metadata that
% gets used in subsequent analyses. I created this so that we can start
% from a clean workspace and do anything after, without having to run
% several sections of a messy analysis script.

load("../datasets/zoopee_pM_OneMode.mat")

% This is included just so when you graph things you can avoid having
% descriptions like adduct names. This is basically the mtabData variable
% but nicer. 
load("nicenames.mat")

% setDefaultFigs
% Load up the color palette of choice. 
load("AlbumMaps.mat", "CP1")

% Notes from lab notebook:
%   Lost 0.5 mL (~1/12) of extract #29 - should not affect L/H ratio
%   During cartridge loading, accidentally put ~2 mL of #22 (ratio not
%   affected) onto #9: check 70 pg/mL standard for contamination from t6
%   ctrl.
%   Same as above but #8 (50 pg/mL) -> #16 (3000 pg/mL)

% Two different sets of metadata, where the first is the sequence list and
% the second is the spreadsheet where Amy and I tracked the experiment. 
Info = '../datasets/mtab_Noah_Zoop2_BC_071423.xlsx';
Info = readtable(Info);
AmyInfo = '../datasets/2023_May_small_boat_metabolomics_filled.xlsx';
AmyInfo = readtable(AmyInfo, 'ReadVariableNames', true);

bigInfo = join(Info(~isnan(Info.bottleNumAmy),:), AmyInfo(~isnan(AmyInfo.BTL_ID),:), 'LeftKeys', {'bottleNumAmy'}, 'RightKeys', {'BTL_ID'});

% Join that data, omitting the multiples for pos+neg modes.
[sI, iBig] = innerjoin(sInfo, bigInfo, ...
    'LeftKeys', 'FileName_neg', 'RightKeys', 'File_Name');

% Reorder metabolite data 
mtabData_pM_Reorder = mtabData_pM(:,iBig);

clear sInfo iBig AmyInfo mtabData_pM bigInfo Info