function [LOD_pM,LOQ_pM,mtabData_pM,MaxStd_pM,mtabElem] = convertMoles(negTransitions, posTransitions, mtabNames, mtabData, LOD, LOQ, MaxStd)

%% Convert from ng added to nM concentrations
% This process involves reloading the transition list files and matching
% names to molecular weights (MW). For DT5, I had to manually correct 
% many, many names. That's why the data file loaded at the beginning 
% contains two variables for mtabNames.
%

posInfo = readtable(posTransitions);
posInfo(posInfo.isParent == 0,:) = [];
MWp = posInfo(:,[1,14:17]); %please confirm column 14 is the MW or the code will not execute properly
negInfo = readtable(negTransitions);
negInfo(negInfo.isParent == 0,:) = [];
MWn = negInfo(:,[1,13:15,17]); %%please confirm column 14 is the MW or the code will not execute properly
MW = [MWp;MWn];
clear MWp MWn
MW = unique(MW, 'rows');
clear posInfo negInfo

% Making both compound name columns into strings and removing the neg/pos
% identifier.
MW.CompoundName = string(MW.CompoundName);
mtabNamesAgnostic = strrep(strrep(mtabNames, ' pos', ''),' neg','');


% Time to index where each unique molecule is found in mtabNames.
[~, iNames] = ismember(mtabNamesAgnostic, MW.CompoundName);
iNames(iNames==0) = [];
% Use those indices to sort MW values.
if sum(mtabNamesAgnostic == MW.CompoundName(iNames)) ~= length(mtabNames)
    disp("name mismatch")
    return
end
MWtoConvert = MW.StdMW(iNames);
mtabElem = MW(iNames, 2:4);
% convert from ng added to ng/L to nM
mtabData_pM = (mtabData.*1000)./MWtoConvert;
LOD_pM = LOD.*1000./MWtoConvert;
LOQ_pM = LOQ.*1000./MWtoConvert;
MaxStd_pM = MaxStd*1000./MWtoConvert;

end