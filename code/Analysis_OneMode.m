% Noah Germolus 25 August 2022
% This is a quick analysis script for the Zoop data. 

clear; clc;  
close all

loadstuff

outdir = '../figs/zoopfigs1';
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

% Information for constructing a color gradient later on. 
colors = [CP1{1};CP1{4};CP1{5};CP1{2}];
npts = 1000;

%% Calculations.
%
%   I want to calculate several things. 
%   -change in concentration from the t0 controls. This will be done by
%       subtracting the mean control at t0 from each concentration, including
%       the additional controls. 
%   -change in inventory: same calculation, using the volume-normalized
%       mtabData_pmol
%   -change at 12h on a dry-mass basis. This will utilize the
%       volume-normalized inventories, but additionally divided by the dry
%       weight of the bugs. 
%

% First: let's actually break the metabolites into elements. 
totC = mtabData_pM_Reorder.*mtabElem.C;
totN = mtabData_pM_Reorder.*mtabElem.N;

% Determine the actual incubation durations. 
sI.duration = sI.TimeStop - sI.TimeStart;

% Create useful indices for controls.
t0ctrl_i = find(sI.Nominal_Duration_h == 0);
t6ctrl_i = find(sI.Nominal_Duration_h == 6 & sI.Species == "CTRL");
t12ctrl_i = find(sI.Nominal_Duration_h == 12 & sI.Species == "CTRL");

% Useful indices for time points and animals.
t0i = sI.Nominal_Duration_h == 0;
t6i = sI.Nominal_Duration_h == 6;
t12i = sI.Nominal_Duration_h==12;
iClio = sI.Species == "C. pyrimidata";
iPx = sI.Species == "PX";
iAmph = sI.Species == "Amphipod (long skinny)";
iEuph = sI.Species == "Euph";
iAnimal = iClio | iPx | iAmph | iEuph;

% Find mean control values. 
meanctrl0 = mean(mtabData_pM_Reorder(:,t0ctrl_i),2);
meanctrl6 = mean(mtabData_pM_Reorder(:,t6ctrl_i),2);
meanctrl12 = mean(mtabData_pM_Reorder(:,t12ctrl_i),2);
stdvctrl0 = std(mtabData_pM_Reorder(:,t0ctrl_i),[],2);
t6ctrl = mtabData_pM_Reorder(:,t6ctrl_i);
stdvctrl6 = std(t6ctrl,[],2);
t12ctrl = mtabData_pM_Reorder(:,t12ctrl_i);
stdvctrl12 = std(t12ctrl,[],2); 
ctrls = [meanctrl0, meanctrl6, meanctrl12];
ste_ctrls = [stdvctrl0, stdvctrl6, stdvctrl12];

% Calculations to determine the inventories and rates to a first degree
% with control subtraction. 
% Subrtract time-matched controls.
mtabData_pM_Reorder_subctrl = mtabData_pM_Reorder;
mtabData_pM_Reorder_subctrl(:,sI.Nominal_Duration_h==0) = mtabData_pM_Reorder_subctrl(:,sI.Nominal_Duration_h==0) - meanctrl0;
mtabData_pM_Reorder_subctrl(:,sI.Nominal_Duration_h==6) = mtabData_pM_Reorder_subctrl(:,sI.Nominal_Duration_h==6) - meanctrl6;
mtabData_pM_Reorder_subctrl(:,sI.Nominal_Duration_h==12) = mtabData_pM_Reorder_subctrl(:,sI.Nominal_Duration_h==12) - meanctrl12;

% Multiply by volume.
mtabData_pmol = (mtabData_pM_Reorder_subctrl' .* (sI.Volume_mL_./1000))';

% Divide by biomass, then...
mtabData_pmol_mgdry = (mtabData_pmol' ./ (sI.dryWeight))';

% ...by time.
mtabData_pmol_mgdry_hr = (mtabData_pmol_mgdry' ./ hours(sI.duration))';

% Times! Yeah, I should probably put some bars on that. 
meantimes = hours([mean(sI.duration(t0i)),...
    mean(sI.duration(t6i)),...
    mean(sI.duration(t12i))]);
rangetimes = hours([range(sI.duration(t0i)),...
    range(sI.duration(t6i)),...
    range(sI.duration(t12i))]);

% Now, the straight concentrations. 
% For this type of data, I'll deal with the P.xiph data first since it has
% the most complicated form. I won't need to repeat the time vector again. 
% This wasn't a good way to go about things--redoing the calcs from a few
% lines above but for the individual species. To-do, maybe: just redo
% subsequent calcs using the larger variables indexed by animal. 

Pxtimes = meantimes(2:end);
Pxtimer = rangetimes(2:end);

Px_pM = [nanmean(mtabData_pM_Reorder(:,t6i&iPx),2),...
    nanmean(mtabData_pM_Reorder(:,t12i&iPx),2)];
Px_stde_pM = [nanstd(mtabData_pM_Reorder(:,t6i&iPx),[],2),...
    nanstd(mtabData_pM_Reorder(:,t12i&iPx),[],2)]; %./sqrt(3);

Amph_pM = nanmean(mtabData_pM_Reorder(:,iAmph),2);
Amph_stde_pM = nanstd(mtabData_pM_Reorder(:,iAmph),[],2);

Clio_pM = nanmean(mtabData_pM_Reorder(:,iClio),2);
Clio_stde_pM = nanstd(mtabData_pM_Reorder(:,iClio),[],2);

Euph_pM = nanmean(mtabData_pM_Reorder(:,iEuph),2);
Euph_stde_pM = nanstd(mtabData_pM_Reorder(:,iEuph),[],2);

% Now, the subtracted versions. 
Px_pM_subctrl = [nanmean(mtabData_pM_Reorder_subctrl(:,t6i&iPx),2),...
    nanmean(mtabData_pM_Reorder_subctrl(:,t12i&iPx),2)];
Px_stde_pM_subctrl = [nanstd(mtabData_pM_Reorder_subctrl(:,t6i&iPx),[],2),...
    nanstd(mtabData_pM_Reorder_subctrl(:,t12i&iPx),[],2)]; %./sqrt(3);

Amph_pM_subctrl = nanmean(mtabData_pM_Reorder_subctrl(:,iAmph),2);
Amph_stde_pM_subctrl = nanstd(mtabData_pM_Reorder_subctrl(:,iAmph),[],2);

Clio_pM_subctrl = nanmean(mtabData_pM_Reorder_subctrl(:,iClio),2);
Clio_stde_pM_subctrl = nanstd(mtabData_pM_Reorder_subctrl(:,iClio),[],2);

Euph_pM_subctrl = nanmean(mtabData_pM_Reorder_subctrl(:,iEuph),2);
Euph_stde_pM_subctrl = nanstd(mtabData_pM_Reorder_subctrl(:,iEuph),[],2);

% Now, the inventory.
Px_pmol = [nanmean(mtabData_pmol(:,t6i&iPx),2),...
    nanmean(mtabData_pmol(:,t12i&iPx),2)];
Px_stde_pmol = [nanstd(mtabData_pmol(:,t6i&iPx),[],2),...
    nanstd(mtabData_pmol(:,t12i&iPx),[],2)]; %./sqrt(3);

Amph_pmol = nanmean(mtabData_pmol(:,iAmph),2);
Amph_stde_pmol = nanstd(mtabData_pmol(:,iAmph),[],2);

Clio_pmol = nanmean(mtabData_pmol(:,iClio),2);
Clio_stde_pmol = nanstd(mtabData_pmol(:,iClio),[],2);

Euph_pmol = nanmean(mtabData_pmol(:,iEuph),2);
Euph_stde_pmol = nanstd(mtabData_pmol(:,iEuph),[],2);

% Now, the time-normalized inventory.
Px_pmol_mgdry_hr = [nanmean(mtabData_pmol_mgdry_hr(:,t6i&iPx),2),...
    nanmean(mtabData_pmol_mgdry_hr(:,t12i&iPx),2)];
Px_stde_pmol_mgdry_hr = [nanstd(mtabData_pmol_mgdry_hr(:,t6i&iPx),[],2),...
    nanstd(mtabData_pmol_mgdry_hr(:,t12i&iPx),[],2)]; %./sqrt(3);

Amph_pmol_mgdry_hr = nanmean(mtabData_pmol_mgdry_hr(:,iAmph),2);
Amph_stde_pmol_mgdry_hr = nanstd(mtabData_pmol_mgdry_hr(:,iAmph),[],2);

Clio_pmol_mgdry_hr = nanmean(mtabData_pmol_mgdry_hr(:,iClio),2);
Clio_stde_pmol_mgdry_hr = nanstd(mtabData_pmol_mgdry_hr(:,iClio),[],2);

Euph_pmol_mgdry_hr = nanmean(mtabData_pmol_mgdry_hr(:,iEuph),2);
Euph_stde_pmol_mgdry_hr = nanstd(mtabData_pmol_mgdry_hr(:,iEuph),[],2);

% Creating a screening process for metabolites that are too high in the
% control or else just bad ones that slipped past the automated QC

iallctrl = sI.Species=="CTRL"; 
mtab_pM_ctrl = mtabData_pM_Reorder(:,iallctrl);
mtab_pM_notctrl = mtabData_pM_Reorder(:,~iallctrl);
mremove = (10.*mean(mtab_pM_ctrl,2,"omitnan")>mean(mtab_pM_notctrl,2,"omitnan"));
rmNames = ["2'deoxyguanosine Na","S-(5'-adenosyl)-L-homocysteine",...
     "adenosine", "serine 2", "4-aminobenzoic acid","cysteine dimer",...
     "ornithine 2"]';
[~, ibad] = ismember(rmNames, mtabNames); mremove(ibad) = 1;

clear rmNames mtab_pM_notctrl mtab_pM_ctrl iallctrl

%% Non-Metric Multidimensional Scaling
% Literally, ignore all the commented stuff. I tested some different
% methods and distance metrics, and more or less landed on Bray-Curtis and
% the construction of a tree with bootstrapping. 

num_boots = 10000;
mD = mtabData_pmol_mgdry_hr(~mremove,:); % Set which normalization to use
% construncting trees based on 
% AnalysisType = "Euclidean";
% AnalysisType = "Jaccard";
% AnalysisType = "BrayCurtis";
%AnalysisType = "BC_BootstrapTree";

[tr,LabelOrder,LabelsOrdered] = bootstrapper(mD,sI,num_boots);

% Heatmap of Metabolites
HeatMapMtabs = mtabData_pmol_mgdry_hr(~mremove,LabelOrder);
HeatMapMtabs(HeatMapMtabs==0) = NaN;
h = heatmap(nicenames(~mremove), LabelsOrdered+ string(1:22)', HeatMapMtabs');
h.ColorLimits = [1,10];
h.ColorScaling = "log";
h.Colormap = flip(cmapper(colors,npts));

%% Eliminating indicators of control contamination and redoing.

[tr_removectrls,LabelOrder2,LabelsOrdered2] = bootstrapper(mD(~mremove,:),sI,num_boots);

% Heatmap of Metabolites
HeatMapMtabs = mtabData_pmol(~mremove,LabelOrder2);
HeatMapMtabs(HeatMapMtabs==0) = NaN;
h2 = heatmap(nicenames(~mremove), LabelsOrdered2+ string(1:22)', HeatMapMtabs');
h2.ColorLimits = [1,10];
h2.ColorScaling = "log";
h2.Colormap = flip(cmapper(colors,npts));

%% Now eliminate controls and do it again. 

[tr_removectrls3,LabelOrder3,LabelsOrdered3] = bootstrapper(mD(~mremove,~iallctrl),sI(~iallctrl,:),num_boots);

% Heatmap of Metabolites
HeatMapMtabs = mtabData_pmol_mgdry_hr(~mremove,~iallctrl);
HeatMapMtabs = HeatMapMtabs(:,LabelOrder3(LabelOrder3~=0));
HeatMapMtabs(HeatMapMtabs==0) = NaN;
h3 = heatmap(nicenames(~mremove), LabelsOrdered3+ string(1:sum(~iallctrl))', HeatMapMtabs');
h3.ColorLimits = [1,10];
h3.ColorScaling = "log";
h3.Colormap = flip(cmapper(colors,npts));

MaxRM = MaxStd_pM(~mremove);
AboveMax = mtabData_pM_Reorder(~mremove,~iallctrl);
AboveMax = AboveMax(:,LabelOrder3(LabelOrder3~=0))>MaxRM;
AboveMax = flip(AboveMax',1);
ax = axes(gcf, "Position", h3.Position, "Units","normalized","Color", "none");
ax.Box = "off";
set(ax, "XTick", [], "YTick", [])
ax.XLim = [0 size(AboveMax,2)]; ax.YLim = [0 size(AboveMax,1)];
[X,Y] = meshgrid(0.5:size(AboveMax,2)-0.5,0.5:size(AboveMax,1)-0.5);
t = text(X(AboveMax), Y(AboveMax),"*");

% switch AnalysisType
%     case "Euclidean"
%         mnan = mtabData_pM_Reorder; % mtabData_pmol_mgdry;
%         mnan(mnan<=0) = NaN;
%         minm = min(min(mnan,[],2,"omitnan"),[],1,"omitnan");
%         c = round(log2(minm)); d = 2^c;
%         mnan(isnan(mnan))=0;
%         ml2 = log2(mnan + d) - c;
%         D = squareform(pdist(ml2',@naneucdist));
%         [Y,eig] = cmdscale(D);
%         pv = eig./sum(eig);
% 
%         % ml2_z = ml2;
%         % ml2_z(isnan(ml2_z)|isinf(ml2_z))=0;
% 
%         Sfc = Y(:,1:2)-mean(Y(:,1:2),1);
%         mlc = ml2'-mean(ml2',1,"omitnan");
%         coefs = zeros(size(mlc,2),2);
%         r2 = zeros(size(mlc,2),1);
%         for ii=1:length(coefs)
%             lm = fitlm(Sfc,mlc(:,ii));
%             coefs(ii,:) = lm.Coefficients.Estimate(2:3)';
%             r2(ii) = lm.Rsquared.Adjusted;
%         end
%         vecs = coefs./sqrt(sum(coefs.^2,2));
%         vecnames = mtabNames(r2>0.7);
%         vecs = vecs(r2>0.7,:);
%         z = zeros(size(vecs,1),1);
% 
%         sI = sI;
%         sI.x = Y(:,1);
%         sI.y = Y(:,2);
% 
%         px = sI.Species=="PX" & sI.Nominal_Duration_h==12;
%         eu = sI.Species=="Euph";
%         am = sI.Species=="Amphipod (long skinny)";
%         cp = sI.Species=="C. pyrimidata";
%         dead = sI.Notes=="DEAD";
% 
% 
%         sc1 = scatter(sI.x(px),sI.y(px),36,CP1{1},"filled");
%         hold on
%         sc2 = scatter(sI.x(am),sI.y(am),36,CP1{2},"filled");
%         sc3 = scatter(sI.x(eu&~dead),sI.y(eu&~dead),36,CP1{3},"filled");
%         sc4 = scatter(sI.x(cp&~dead),sI.y(cp&~dead),36,CP1{4},"filled");
%         sc6 = scatter(sI.x(cp&dead),sI.y(cp&dead),36,CP1{4}, "HandleVisibility","off");
%         sc5 = scatter(sI.x(eu&dead),sI.y(eu&dead),36, CP1{3}, "HandleVisibility","off");
%         q = quiver(z,z,100*vecs(:,1),100*vecs(:,2),"HandleVisibility","off");
%         t = text(100*vecs(:,1),100*vecs(:,2),vecnames,"HandleVisibility","off");
%         t2 = text(sI.x(cp&dead)+1,sI.y(cp&dead),"dead", "HandleVisibility","off");
%         t3 = text(sI.x(eu&dead)+1,sI.y(eu&dead),"dead", "HandleVisibility","off");
%         legend({'Copepod','Amphipod','Euphausiid','Pteropod'},"Location","best")
%         xlabel(["PC1 ("+string(round(100*pv(1),2))+"%)"]);
%         ylabel(["PC2 ("+string(round(100*pv(2),2))+"%)"]);
%         title("PCA of log_2 Normalized Metabolite Data")
% 
%         % Dendrogram of the Euclidean Distance
% 
%         Labels = [string(sI.Species) + " " + string(round(hours(sI.duration))) + "h " + string(sI.Notes)];
%         Z=linkage(D,'average');
%         [H,T,outperm] = dendrogram(Z,size(sI,1),'Labels',Labels,'orientation','left');
%         LabelOrder = flip(outperm);
% 
%     case "Jaccard"
%         % Jaccard Distance with many random initial configs
%         numTries = 100;
%         mn = mtabData_pM_Reorder; mn(isnan(mn))=0;
%         D = squareform(pdist(mn', 'jaccard'));
%         opts = statset('MaxIter',10000,'Display','off','TolX',1e-5);
%         Y0 = zeros(size(mn,2),2,numTries);
%         stress0 = zeros(numTries,1);
%         for i=1:numTries
%             try
%                 [Y0(:,:,ii),stress0(ii)] = mdscale(D, 2, 'Start','random','Options',opts);
%             catch error
%                 Y0(:,:,ii) = NaN;
%                 stress0(ii) = NaN;
%             end
% 
%         end
%         Y = 1e4.*mean(Y0,3, 'omitnan');
%         stress = mean(stress0, 1, 'omitnan');
% 
%         % ml2_z = ml2;
%         % ml2_z(isnan(ml2_z)|isinf(ml2_z))=0;
% 
%         sI = sI;
%         sI.x = Y(:,1);
%         sI.y = Y(:,2);
% 
%         px = sI.Species=="PX" & sI.Nominal_Duration_h==12;
%         eu = sI.Species=="Euph";
%         am = sI.Species=="Amphipod (long skinny)";
%         cp = sI.Species=="C. pyrimidata";
%         dead = sI.Notes=="DEAD";
% 
%         mnr = mn(sum(mn,2)>mean(sum(mn,2)),:);
%         nr = mtabNames(sum(mn,2)>mean(sum(mn,2)));
%         Sfc = Y(:,1:2)-mean(Y(:,1:2),1);
%         mlc = [mnr',sI.dryWeight];
%         mlc = mlc-mean(mlc,1,"omitnan");
%         coefs = zeros(size(mlc,2),2);
%         r2 = zeros(size(mlc,2),1);
%         for ii=1:size(coefs,1)
%             lm = fitlm(Sfc,mlc(:,ii));
%             coefs(ii,:) = lm.Coefficients.Estimate(2:3)';
%             r2(ii) = lm.Rsquared.Adjusted;
%         end
%         vecs = r2.*coefs./sqrt(sum(coefs.^2,2));
%         vecnames = [nr;"Dry Weight"]; %mtabNames(r2>0.7);
%         %vecs = vecs(r2>0.7,:);
%         z = zeros(size(vecs,1),1);
% 
% 
%         sc1 = scatter(sI.x(px),sI.y(px),36,CP1{1},"filled");
%         hold on
%         sc2 = scatter(sI.x(am),sI.y(am),36,CP1{2},"filled");
%         sc3 = scatter(sI.x(eu&~dead),sI.y(eu&~dead),36,CP1{3},"filled");
%         sc4 = scatter(sI.x(cp&~dead),sI.y(cp&~dead),36,CP1{4},"filled");
%         sc6 = scatter(sI.x(cp&dead),sI.y(cp&dead),36,CP1{4}, "HandleVisibility","off");
%         sc5 = scatter(sI.x(eu&dead),sI.y(eu&dead),36, CP1{3}, "HandleVisibility","off");
%         sc7 = scatter(sI.x(t12ctrl_i), sI.y(t12ctrl_i), 36, CP1{5}, 'filled');
%         q = quiver(z,z,vecs(:,1),vecs(:,2),"HandleVisibility","off");
%         t = text(vecs(:,1),vecs(:,2),vecnames,"HandleVisibility","off");
%         t2 = text(sI.x(cp&dead)+0.1,sI.y(cp&dead),"dead", "HandleVisibility","off");
%         t3 = text(sI.x(eu&dead)+0.1,sI.y(eu&dead),"dead", "HandleVisibility","off");
%         legend({'Copepod','Amphipod','Euphausiid','Pteropod','control'},"Location","best")
%         xlabel("NMDS1");
%         ylabel("NMDS2");
%         title("NMDS of Metabolite Data, Stress = " + string(stress))
% 
%         % Dendrogram for Jaccard
% 
%         Labels = [string(sI.Species) + " " + string(round(hours(sI.duration))) + "h " + string(sI.Notes)];
%         Z=linkage(D,'average');
%         [H,T,outperm] = dendrogram(Z,size(sI,1),'Labels',Labels,'orientation','left');
%         LabelOrder = flip(outperm);
% 
%     case "BrayCurtis"
% 
        % %% Bray-Curtis with multiple random starts
        % numTries = 10;
        % mn = mD(~mremove,:); mn(isnan(mn))=0; %(~mremove,~iallctrl)
        % D = f_braycurtis(mn);
        % opts = statset('MaxIter',10000,'Display','off','TolX',1e-5);
        % Y0 = zeros(size(mn,2),2,numTries);
        % stress0 = zeros(numTries,1);
        % for ii=1:numTries
        %     [Y0(:,:,ii),stress0(ii)] = mdscale(D, 2, 'Start','random','Options',opts);
        % end
        % Y = 1e4.*mean(Y0,3);
        % stress = mean(stress0);
        % 
        % % ml2_z = ml2;
        % % ml2_z(isnan(ml2_z)|isinf(ml2_z))=0;
        % 
        % sI = sI;
        % % sI.x(~iallctrl) = Y(:,1);
        % % sI.y(~iallctrl) = Y(:,2);
        % sI.x = Y(:,1);
        % sI.y = Y(:,2);
        % 
        % px = sI.Species=="PX" & sI.Nominal_Duration_h==12;
        % eu = sI.Species=="Euph";
        % am = sI.Species=="Amphipod (long skinny)";
        % cp = sI.Species=="C. pyrimidata";
        % dead = sI.Notes=="DEAD";
        % 
        % mnr = mn(sum(mn,2)>mean(sum(mn,2)),:);
        % nr = mtabNames(~mremove,:);
        % nr = nr(sum(mn,2)>mean(sum(mn,2)));
        % Sfc = Y(:,1:2)-mean(Y(:,1:2),1);
        % mlc = [mnr',sI.dryWeight]; %(~iallctrl)
        % mlc = mlc-mean(mlc,1,"omitnan");
        % coefs = zeros(size(mlc,2),2);
        % r2 = zeros(size(mlc,2),1);
        % for ii=1:size(coefs,1)
        %     lm = fitlm(Sfc,mlc(:,ii));
        %     coefs(ii,:) = lm.Coefficients.Estimate(2:3)';
        %     r2(ii) = lm.Rsquared.Adjusted;
        % end
        % vecs = 1000.*r2.*coefs./sqrt(sum(coefs.^2,2));
        % vecnames = [nr;"Dry Weight"]; %mtabNames(r2>0.7);
        % [~, InterestIndex] = ismember(["Dry Weight", "arginine", "lysine 2", "serine","putrescine"],vecnames);
        % if sum(InterestIndex==0)>0
        %     InterestIndex(:,InterestIndex==0) = [];
        % end
        % %vecs = vecs(r2>0.7,:);
        % z = zeros(size(vecs,1),1);
        % 
        % sc1 = scatter(sI.x(px),sI.y(px),100,CP1{1},"filled"); %(~iallctrl & px)
        % hold on
        % sc2 = scatter(sI.x(am),sI.y(am),100,CP1{2},"filled");
        % sc3 = scatter(sI.x(eu&~dead),sI.y(eu&~dead),100,CP1{3},"filled");
        % sc4 = scatter(sI.x(cp&~dead),sI.y(cp&~dead),100,CP1{4},"filled");
        % sc6 = scatter(sI.x(cp&dead),sI.y(cp&dead),100,CP1{4}, "HandleVisibility","off");
        % sc5 = scatter(sI.x(eu&dead),sI.y(eu&dead),100, CP1{3}, "HandleVisibility","off");
        % sc7 = scatter(sI.x(t12ctrl_i), sI.y(t12ctrl_i), 100, CP1{5}, 'filled');
        % q = quiver(z,z,vecs(:,1),vecs(:,2),"HandleVisibility","off", "LineWidth",2,"Color", "k", "AutoScale", "off");
        % t = text(1.1.*vecs(InterestIndex,1),1.1.*vecs(InterestIndex,2),vecnames(InterestIndex),"HandleVisibility","off");
        % %t2 = text(sI.x(cp&dead)+0.1,sI.y(cp&dead),"dead", "HandleVisibility","off");
        % %t3 = text(sI.x(eu&dead)+0.1,sI.y(eu&dead),"dead", "HandleVisibility","off");
        % legend({'Copepod','Amphipod','Euphausiid','Pteropod','Control'},"Location","best")%,'control'},"Location","best")
        % xlabel("NMDS1");
        % ylabel("NMDS2");
        % % ylim([-250 1000])
        % title("NMDS of Metabolite Data, Bray-Curtis Dissimilarity, Stress = " + string(stress))
% 
%         % Dendrogram for Bray-Curtis
% 
%         %distfunc = @(x,dnu) squareform(f_braycurtis(x,dnu));
%         Labels = [string(sI.Species(~iallctrl)) + " " + string(round(hours(sI.duration(~iallctrl)))) + "h " + string(sI.Notes(~iallctrl))];
%         Z=linkage(D,'average');
%         [H,T,outperm] = dendrogram(Z,size(sI(~iallctrl,:),1),'Labels',Labels,'orientation','left');
%         LabelOrder = flip(outperm);
% % 
%     case "BC_BootstrapTree"
%         num_boots = 1000
%         [tr,LabelOrder,LabelsOrdered] = bootstrapper(mtabData_pM_Reorder,sI,num_boots);
% end

%% Heatmap of just the animals. 

HeatMapMtabs = mtabData_pmol_mgdry_hr(~mremove, iAnimal);
relevantInfo = sI(iAnimal, :);
Labels = string(relevantInfo.Species) + " " + string(round(hours(relevantInfo.duration))) +...
    "h " + string(relevantInfo.Notes) + " " + string(1:size(relevantInfo,1))';
[OrderedLabels,ia] = sort(Labels);
HeatMapMtabs = HeatMapMtabs(:,ia);
HeatMapMtabs(HeatMapMtabs==0) = NaN;
h3 = heatmap(nicenames(~mremove), OrderedLabels + string(1:sum(~iallctrl))', HeatMapMtabs');
h3.ColorLimits = [1,10];
h3.ColorScaling = "log";
h3.Colormap = flip(cmapper(colors,npts));
MaxRM = MaxStd_pM(~mremove);
AboveMax = mtabData_pM_Reorder(~mremove,~iallctrl);
AboveMax = AboveMax(:,ia)>MaxRM;
AboveMax = flip(AboveMax',1);
ax = axes(gcf, "Position", h3.Position, "Units","normalized","Color", "none");
ax.Box = "off";
set(ax, "XTick", [], "YTick", [])
ax.XLim = [0 size(AboveMax,2)]; ax.YLim = [0 size(AboveMax,1)];
[X,Y] = meshgrid(0.5:size(AboveMax,2)-0.5,0.5:size(AboveMax,1)-0.5);
t = text(X(AboveMax), Y(AboveMax),"*");

%% Boxplot

boxplotmtabs = (mtabData_pM_Reorder(~mremove,:)' .* (sI.Volume_mL_./1000))';


figure 
b = boxplot(boxplotmtabs,sI.Species, "PlotStyle", "compact",...
    "Colors", [CP1{3}; CP1{1}; CP1{2}; CP1{5}; CP1{4}]);
ax = gca; ax.YScale = "log";
ax.XTick = [1,2,3,4,5];
ax.XTickLabel = {"C. pyrimidata", "Control","P. xiphias", "Euphausiid", "Amphipod"};
ax.YLabel.String = "pmol mtab at 12 h";
title("All Metabolites by Sample Type")

%% Scatter

scattertabs = mtabData_pmol(~mremove,t12i)';
scattertabs(isnan(scattertabs))= 0;
rI = sI(sI.Nominal_Duration_h==12,:);

med = median(scattertabs,2,"includemissing");
[~,ia] = sort(med);
scattertabs = scattertabs(ia,:);
scatternames = nicenames(~mremove); scatternames = scatternames(ia);

xCt = repmat(0.6:1:(0.6+length(nicenames(~mremove))-1),1,length(rI.Species(rI.Species=="CTRL")))';
xPx = repmat(0.8:1:(0.8+length(nicenames(~mremove))-1),1,length(rI.Species(rI.Species=="PX")))';
xCp = repmat(1:1:(1+length(nicenames(~mremove))-1),1,length(rI.Species(rI.Species=="C. pyrimidata")))';
xAm = repmat(1.2:1:(1.2+length(nicenames(~mremove))-1),1,length(rI.Species(rI.Species=="Amphipod (long skinny)")))';
xEu = repmat(1.4:1:(1.4+length(nicenames(~mremove))-1),1,length(rI.Species(rI.Species=="Euph")))';
xL = 1.5:1:(1.5+length(nicenames(~mremove))-2);

bmmCt = reshape(scattertabs(:,rI.Species == "CTRL"),length(xCt),1);
bmmPx = reshape(scattertabs(:,rI.Species == "PX"),length(xPx),1);
bmmCp = reshape(scattertabs(:,rI.Species == "C. pyrimidata"),length(xCp),1);
bmmAm = reshape(scattertabs(:,rI.Species == "Amphipod (long skinny)"),length(xAm),1);
bmmEu = reshape(scattertabs(:,rI.Species == "Euph"),length(xEu),1);

figure
xline(xL, "LineStyle","--", "LineWidth",0.5, "Color",[.5,.5,.5],"HandleVisibility","off")
hold on
sc1 = scatter(xCt,bmmCt, 40, CP1{1},"o");
sc2 = scatter(xPx,bmmPx, 40, CP1{2},'filled',"^");
sc3 = scatter(xCp,bmmCp, 40, CP1{3},'filled',"square");
sc4 = scatter(xAm,bmmAm, 40, CP1{4},'filled',"diamond");
sc5 = scatter(xEu,bmmEu, 40, CP1{5},'filled',"hexagram");
ax = gca;
ax.YLabel.String = "metabolite present at 12 h, pmol";
set(ax,"XLim",[0 1.5+length(scatternames)], "YScale","log");
set(ax,"XTick",1:length(scatternames), "XTickLabels", scatternames,"XTickLabelRotationMode","auto")
set(ax, "Box", "off")
legend({"Control", "{\it P. xiphias}", "{\it C. Pyrimidata}", "Amphipod", "Euphausiid"})
title("Control-Screened Metabolites, Volume-Normalized, Ordered by Median")

%% The same but with an alternate orientation
figure
yline(xL, "LineStyle","--", "LineWidth",0.5, "Color",[.5,.5,.5])
hold on
sc1 = scatter(bmmCt,xCt, 40, CP1{1},"o");
sc2 = scatter(bmmPx,xPx, 40, CP1{2},'filled',"^");
sc3 = scatter(bmmCp,xCp, 40, CP1{3},'filled',"square");
sc4 = scatter(bmmAm,xAm, 40, CP1{4},'filled',"diamond");
sc5 = scatter(bmmEu,xEu, 40, CP1{5},'filled',"hexagram");
ax = gca;
ax.XLabel.String = "metabolite present at 12 h, pmol";
set(ax,"YLim",[0 1.5+length(scatternames)], "XScale","log");
set(ax,"YTick",1:length(scatternames), "YTickLabels", scatternames,"YTickLabelRotationMode","auto")
set(ax, "Box", "off", "XAxisLocation", "top")

%% Now a scatter with controls subtracted. 

scattertabs = mtabData_pmol(~mremove,t12i)';
scattertabs(isnan(scattertabs))= 0;
rI = sI(sI.Nominal_Duration_h==12,:);
ictrl = rI.Species=="CTRL";
rI = rI(rI.Species~="CTRL",:);
meanctrlsc = mean(scattertabs(:,ictrl),2,"omitnan");
scattertabs = scattertabs(:,~ictrl)-meanctrlsc;
scattertabs(scattertabs<0)=0;

med = median(scattertabs,2,"includemissing");
[~,ia] = sort(med);
scattertabs = scattertabs(ia,:);
scatternames = nicenames(~mremove); scatternames = scatternames(ia);

%xCt = repmat(0.6:1:(0.6+length(nicenames(~mremove))-1),1,length(rI.Species(rI.Species=="CTRL")))';
xPx = repmat(0.8:1:(0.8+length(nicenames(~mremove))-1),1,length(rI.Species(rI.Species=="PX")))';
xCp = repmat(1:1:(1+length(nicenames(~mremove))-1),1,length(rI.Species(rI.Species=="C. pyrimidata")))';
xAm = repmat(1.2:1:(1.2+length(nicenames(~mremove))-1),1,length(rI.Species(rI.Species=="Amphipod (long skinny)")))';
xEu = repmat(1.4:1:(1.4+length(nicenames(~mremove))-1),1,length(rI.Species(rI.Species=="Euph")))';
xL = 1.5:1:(1.5+length(nicenames(~mremove))-2);

%bmmCt = reshape(scattertabs(:,rI.Species == "CTRL"),length(xCt),1);
bmmPx = reshape(scattertabs(:,rI.Species == "PX"),length(xPx),1);
bmmCp = reshape(scattertabs(:,rI.Species == "C. pyrimidata"),length(xCp),1);
bmmAm = reshape(scattertabs(:,rI.Species == "Amphipod (long skinny)"),length(xAm),1);
bmmEu = reshape(scattertabs(:,rI.Species == "Euph"),length(xEu),1);

figure
xline(xL, "LineStyle","--", "LineWidth",0.5, "Color",[.5,.5,.5],"HandleVisibility","off")
hold on
%sc1 = scatter(xCt,bmmCt, 40, CP1{1},"o");
sc2 = scatter(xPx,bmmPx, 40, CP1{2},'filled',"^");
sc3 = scatter(xCp,bmmCp, 40, CP1{3},'filled',"square");
sc4 = scatter(xAm,bmmAm, 40, CP1{4},'filled',"diamond");
sc5 = scatter(xEu,bmmEu, 40, CP1{5},'filled',"hexagram");
ax = gca;
ax.YLabel.String = "metabolite present at 12 h, pmol";
set(ax,"XLim",[0 1.5+length(scatternames)], "YScale","log");
set(ax,"XTick",1:length(scatternames), "XTickLabels", scatternames,"XTickLabelRotationMode","auto")
set(ax, "Box", "off")
legend({"{\it P. xiphias}", "{\it C. pyrimidata}", "Amphipod", "Euphausiid"})
title("Control-{\it Subtracted} Metabolites, Volume-Normalized, Ordered by Median")

%% Now a scatter with controls subtracted and mass-normalized, dead removed. 

scattertabs = mtabData_pmol(~mremove,t12i)';
scattertabs(isnan(scattertabs))= 0;
rI = sI(sI.Nominal_Duration_h==12,:);
ictrl = rI.Species=="CTRL";
rI = rI(rI.Species~="CTRL",:);
meanctrlsc = mean(scattertabs(:,ictrl),2,"omitnan");
scattertabs = scattertabs(:,~ictrl)-meanctrlsc;
scattertabs(scattertabs<0)=0;
scattertabs = scattertabs./rI.dryWeight'./hours(rI.duration)';
deadCp = rI.Notes(rI.Species=="C. pyrimidata")=="DEAD";
deadEu = rI.Notes(rI.Species=="Euph")=="DEAD";
deadCp = repelem(deadCp',1,length(nicenames(~mremove)));
deadEu = repelem(deadEu',1,length(nicenames(~mremove)));

med = median(scattertabs,2,"includemissing");
[~,ia] = sort(med);
scattertabs = scattertabs(ia,:);
scatternames = nicenames(~mremove); scatternames = scatternames(ia);

xPx = repmat(0.8:1:(0.8+length(nicenames(~mremove))-1),1,length(rI.Species(rI.Species=="PX")))';
xCp = repmat(1:1:(1+length(nicenames(~mremove))-1),1,length(rI.Species(rI.Species=="C. pyrimidata")))';
xAm = repmat(1.2:1:(1.2+length(nicenames(~mremove))-1),1,length(rI.Species(rI.Species=="Amphipod (long skinny)")))';
xEu = repmat(1.4:1:(1.4+length(nicenames(~mremove))-1),1,length(rI.Species(rI.Species=="Euph")))';
xL = 1.5:1:(1.5+length(nicenames(~mremove))-2);

bmmPx = reshape(scattertabs(:,rI.Species == "PX"),length(xPx),1);
bmmCp = reshape(scattertabs(:,rI.Species == "C. pyrimidata"),length(xCp),1);
bmmAm = reshape(scattertabs(:,rI.Species == "Amphipod (long skinny)"),length(xAm),1);
bmmEu = reshape(scattertabs(:,rI.Species == "Euph"),length(xEu),1);

figure
xline(xL, "LineStyle","--", "LineWidth",0.5, "Color",[.5,.5,.5],"HandleVisibility","off")
hold on
sc2 = scatter(xPx,bmmPx, 40, CP1{2},'filled',"^");
sc3 = scatter(xCp(~deadCp),bmmCp(~deadCp), 40, CP1{3},'filled',"square");
sc4 = scatter(xAm,bmmAm, 40, CP1{4},'filled',"diamond");
sc5 = scatter(xEu(~deadEu),bmmEu(~deadEu), 40, CP1{5},'filled',"hexagram");
sc6 = scatter(xCp(deadCp),bmmCp(deadCp), 40, CP1{3},"square");
sc7 = scatter(xEu(deadEu),bmmEu(deadEu), 40, CP1{5},"hexagram");
ax = gca;
ax.YLabel.String = "pmol h^{-1} mg^{-1}";
set(ax,"XLim",[0 1.5+length(scatternames)], "YScale","log");
set(ax,"XTick",1:length(scatternames), "XTickLabels", scatternames,"XTickLabelRotationMode","auto")
set(ax, "Box", "off")
legend({"{\it P. xiphias}", "{\it C. pyrimidata}", "Amphipod", "Euphausiid"})
title("Metabolites, Biomass-Normalized, Ordered by Median")

%% Now a scatter with controls subtracted and mass-normalized, dead removed. 
% AND ELIMINATE ANYTHING UNDER 2.1667

scattertabs = mtabData_pmol(~mremove,t12i)';
scattertabs(isnan(scattertabs))= 0;
rI = sI(sI.Nominal_Duration_h==12,:);
ictrl = rI.Species=="CTRL";
rI = rI(rI.Species~="CTRL",:);
meanctrlsc = mean(scattertabs(:,ictrl),2,"omitnan");
scattertabs = scattertabs(:,~ictrl)-meanctrlsc;
scattertabs(scattertabs<0)=0;
scattertabs = scattertabs./rI.dryWeight'./hours(rI.duration)';
%scattertabs(scattertabs<2.1667) = 0;
below = sum(scattertabs<2.1667,2)==0;
deadCp = rI.Notes(rI.Species=="C. pyrimidata")=="DEAD";
deadEu = rI.Notes(rI.Species=="Euph")=="DEAD";
deadCp = repelem(deadCp',1,sum(~below));
deadEu = repelem(deadEu',1,sum(~below));

med = median(scattertabs,2,"includemissing");
[OrderedMedians,ia] = sort(med);
scattertabs = scattertabs(ia(~below),:);
scatternames = nicenames(~mremove); 
scatternames = scatternames(ia(~below));

xPx = repmat(0.8:1:(0.8+sum(~below)-1),1,length(rI.Species(rI.Species=="PX")))';
xCp = repmat(1:1:(1+sum(~below)-1),1,length(rI.Species(rI.Species=="C. pyrimidata")))';
xAm = repmat(1.2:1:(1.2+sum(~below)-1),1,length(rI.Species(rI.Species=="Amphipod (long skinny)")))';
xEu = repmat(1.4:1:(1.4+sum(~below)-1),1,length(rI.Species(rI.Species=="Euph")))';
xL = 1.5:1:(1.5+sum(~below)-2);

bmmPx = reshape(scattertabs(:,rI.Species == "PX"),length(xPx),1);
bmmCp = reshape(scattertabs(:,rI.Species == "C. pyrimidata"),length(xCp),1);
bmmAm = reshape(scattertabs(:,rI.Species == "Amphipod (long skinny)"),length(xAm),1);
bmmEu = reshape(scattertabs(:,rI.Species == "Euph"),length(xEu),1);

figure
xline(xL, "LineStyle","--", "LineWidth",0.5, "Color",[.5,.5,.5],"HandleVisibility","off")
hold on
sc2 = scatter(xPx,bmmPx, 40, CP1{2},'filled',"^");
sc3 = scatter(xCp(~deadCp),bmmCp(~deadCp), 40, CP1{3},'filled',"square");
sc4 = scatter(xAm,bmmAm, 40, CP1{4},'filled',"diamond");
sc5 = scatter(xEu(~deadEu),bmmEu(~deadEu), 40, CP1{5},'filled',"hexagram");
sc6 = scatter(xCp(deadCp),bmmCp(deadCp), 40, CP1{3},"square");
sc7 = scatter(xEu(deadEu),bmmEu(deadEu), 40, CP1{5},"hexagram");
ax = gca;
ax.YLabel.String = "pmol h^{-1} mg^{-1}";
set(ax,"XLim",[0 1.5+length(scatternames)], "YScale","log");
set(ax,"XTick",1:length(scatternames), "XTickLabels", scatternames,"XTickLabelRotationMode","auto")
set(ax, "Box", "off", "YLim", [2.1667, 1e4])
legend({"{\it P. xiphias}", "{\it C. pyrimidata}", "Amphipod", "Euphausiid"})
title("Metabolites, Biomass-Normalized, Ordered by Median")

%% The same as above but medians by species. 

scattertabs = mtabData_pmol(~mremove,t12i)';
scattertabs(isnan(scattertabs))= 0;
rI = sI(sI.Nominal_Duration_h==12,:);
ictrl = rI.Species=="CTRL";
rI = rI(~ictrl,:);
scattertabs = scattertabs(~ictrl,:);
scattertabs(scattertabs<0)=0;
% ia12 = rI.Nominal_Duration_h==12&(rI.Species=="PX"|rI.Species=="C. pyrimidata"|rI.Species=="Euph"|rI.Species=="Amphipod (long skinny)");
scattertabs = scattertabs./rI.dryWeight./hours(rI.duration);
% clear ia12

subplot(2,2,1)
% P. xiphias plot.
pxtabs = scattertabs(rI.Species == "PX",:);
below = sum(pxtabs>2.1667,1)==0;
med = median(pxtabs,1,"includemissing");
[~,ia] = sort(med);
pxtabs = pxtabs(:,ia(~below));
pxnames = nicenames(~mremove); 
pxnames = pxnames(ia(~below));
xPx = repmat(1:sum(~below),1,length(rI.Species(rI.Species=="PX")))';
bmmPx = reshape(pxtabs,length(xPx),1);
sc1 = scatter(xPx,bmmPx, 40, CP1{2},'filled',"^");
ax = gca;
ax.YLabel.String = "pmol h^{-1} mg^{-1}";
set(ax,"XLim",[0 1+length(pxnames)], "YScale","log");
set(ax,"XTick",1:length(pxnames), "XTickLabels", pxnames,"XTickLabelRotation",45)
hold on
xL = 1.5:1:(1.5+sum(~below)-2);
xline(xL, "LineStyle","--", "LineWidth",0.5, "Color",[.5,.5,.5],"HandleVisibility","off")
title("Pleuromamma xiphias", "FontAngle","italic")

subplot(2,2,2)
% C. py plot.
cptabs = scattertabs(rI.Species == "C. pyrimidata",:);
below = sum(cptabs>2.1667,1)==0;
med = median(cptabs,1,"includemissing");
[~,ia] = sort(med);
tempI = rI(rI.Species=="C. pyrimidata",:);
idead = tempI.Notes=="DEAD";
idead = repmat(idead',1,sum(~below))';
cptabs = cptabs(:,ia(~below));
cpnames = nicenames(~mremove); 
cpnames = cpnames(ia(~below));
xCp = repelem(1:sum(~below),1,length(rI.Species(rI.Species=="C. pyrimidata")))';
bmmCp = reshape(cptabs,length(xCp),1);
sc2 = scatter(xCp(~idead),bmmCp(~idead), 40, CP1{3},'filled',"square");
hold on
sc3 = scatter(xCp(idead),bmmCp(idead), 40, CP1{3},"square");
ax = gca;
ax.YLabel.String = "pmol h^{-1} mg^{-1}";
set(ax,"XLim",[0 1+length(cpnames)], "YScale","log");
set(ax,"XTick",1:length(cpnames), "XTickLabels", cpnames,"XTickLabelRotation",45)
hold on
xL = 1.5:1:(1.5+sum(~below)-2);
xline(xL, "LineStyle","--", "LineWidth",0.5, "Color",[.5,.5,.5],"HandleVisibility","off")
clear tempI idead
title("Clio pyrimidata", "FontAngle","italic")

subplot(2,2,3)
% Euph plot.
eutabs = scattertabs(rI.Species == "Euph",:);
below = sum(eutabs>2.1667,1)==0;
med = median(eutabs,1,"includemissing");
[~,ia] = sort(med);
tempI = rI(rI.Species=="Euph",:);
idead = tempI.Notes=="DEAD";
idead = repmat(idead',1,sum(~below))';
eutabs = eutabs(:,ia(~below));
eunames = nicenames(~mremove); 
eunames = eunames(ia(~below));
xEu = repelem(1:sum(~below),1,length(rI.Species(rI.Species=="Euph")))';
bmmEu = reshape(eutabs,length(xEu),1);
sc2 = scatter(xEu(~idead),bmmEu(~idead), 40, CP1{5},'filled',"hexagram");
hold on
sc3 = scatter(xEu(idead),bmmEu(idead), 40, CP1{5},"hexagram");
ax = gca;
ax.YLabel.String = "pmol h^{-1} mg^{-1}";
set(ax,"XLim",[0 1+length(eunames)], "YScale","log");
set(ax,"XTick",1:length(eunames), "XTickLabels", eunames,"XTickLabelRotation",45)
hold on
xL = 1.5:1:(1.5+sum(~below)-2);
xline(xL, "LineStyle","--", "LineWidth",0.5, "Color",[.5,.5,.5],"HandleVisibility","off")
clear tempI idead
title("Euphausiid")

subplot(2,2,4)
% Amph plot
amtabs = scattertabs(rI.Species == "Amphipod (long skinny)",:);
below = sum(amtabs>2.1667,1)==0;
med = median(amtabs,1,"includemissing");
[~,ia] = sort(med);
amtabs = amtabs(:,ia(~below));
amnames = nicenames(~mremove); 
amnames = amnames(ia(~below));
xAm = repelem(1:sum(~below),1,length(rI.Species(rI.Species=="Amphipod (long skinny)")))';
bmmAm = reshape(amtabs,length(xAm),1);
sc2 = scatter(xAm,bmmAm, 40, CP1{4},'filled',"diamond");
hold on
ax = gca;
ax.YLabel.String = "pmol h^{-1} mg^{-1}";
set(ax,"XLim",[0 1+length(amnames)], "YScale","log");
set(ax,"XTick",1:length(amnames), "XTickLabels", amnames,"XTickLabelRotation",45)
hold on
xL = 1.5:1:(1.5+sum(~below)-2);
xline(xL, "LineStyle","--", "LineWidth",0.5, "Color",[.5,.5,.5],"HandleVisibility","off")
clear tempI idead
title("Amphipod")

%% Now, account for instraspecies vs. total variance.

scattertabs = mtabData_pmol(~mremove,t12i)';
scattertabs(isnan(scattertabs))= 0;
rI = sI(sI.Nominal_Duration_h==12,:);
ictrl = rI.Species=="CTRL";
rI = rI(~ictrl,:);
scattertabs = scattertabs(~ictrl,:);
scattertabs(scattertabs<0)=0;
scattertabs = scattertabs./rI.dryWeight./hours(rI.duration);

iClio = rI.Species == "C. pyrimidata";
iPx = rI.Species == "PX";
iAmph = rI.Species == "Amphipod (long skinny)";
iEuph = rI.Species == "Euph";

vC = std(scattertabs(iClio,:),[],1);
vP = std(scattertabs(iPx,:),[],1);
vA = std(scattertabs(iAmph,:),[],1);
vE = std(scattertabs(iEuph,:),[],1);

vT = std(scattertabs,[],1);

vCT = vC./vT; vPT = vP./vT; vAT = vA./vT; vET = vE./vT;
vCT(vCT==0)=NaN;
vPT(vPT==0)=NaN;
vET(vET==0)=NaN;
vAT(vAT==0)=NaN;

hLim = 0.2;
vTT = [vCT;vPT;vET;vAT];
allnan = sum(isnan(vTT),1)==size(vTT,1);
aboveLim = sum(vTT>hLim,1)==size(vTT,1);

scatternames = nicenames(~mremove);
xL = 1.5:1:(1.5+sum(~mremove)-2);
xPx = 0.8:1:(0.8+sum(~mremove)-1);
xCp = 1:1:(1+sum(~mremove)-1);
xAm = 1.2:1:(1.2+sum(~mremove)-1);
xEu = 1.4:1:(1.4+sum(~mremove)-1);

figure
xline(xL(1:end-sum(aboveLim)), "LineStyle","-", "LineWidth",0.5, "Color",[.2,.2,.2],"HandleVisibility","off")
hold on
sc2 = scatter(xPx(1:end-sum(aboveLim)),vPT(~aboveLim), 40, CP1{2},'filled',"^");
sc4 = scatter(xAm(1:end-sum(aboveLim)),vAT(~aboveLim), 40, CP1{4},'filled',"diamond");
sc6 = scatter(xCp(1:end-sum(aboveLim)),vCT(~aboveLim), 40, CP1{3},"filled","square");
sc7 = scatter(xEu(1:end-sum(aboveLim)),vET(~aboveLim), 40, CP1{5},"filled","hexagram");
ax = gca;
ax.YLabel.String = "\sigma_{species} \sigma_{all}^{-1}";
set(ax,"XLim",[0 1.5+length(scatternames(~aboveLim))], "YLim", [0 hLim]);
set(ax,"XTick",1:length(scatternames(~aboveLim)), "XTickLabels", scatternames(~aboveLim),"XTickLabelRotationMode","auto")
set(ax, "Box", "off")
legend({"{\it P. xiphias}", "{\it C. pyrimidata}", "Amphipod", "Euphausiid"})
title("Excreted metabolite variance as a fraction of total variance. (<20%)")

%% Finally, an accounting of elements. 

load("DOC_DON.mat")

elmC = totC(:,sI.Nominal_Duration_h==12);
elmN = totN(:,sI.Nominal_Duration_h==12);
elmC(isnan(elmC))= 0; elmN(isnan(elmN))= 0;
rI = sI(sI.Nominal_Duration_h==12,:);

ictrldead = (rI.Species=="CTRL" | rI.Notes=="DEAD");
ictrl = rI.Species=="CTRL";
rI = rI(~ictrldead,:);

meanctrlsc = mean(elmC(:,ictrl),2,"omitnan");
meanctrlsn = mean(elmN(:,ictrl),2,"omitnan");
elmC = elmC(:,~ictrldead)-meanctrlsc;
elmN = elmN(:,~ictrldead)-meanctrlsn;
elmC(elmC<0)=0;elmN(elmN<0)=0;
itaurine = (mtabNames=="taurine"); iglycine = (mtabNames=="glycine");
Ctau = elmC(itaurine,:); Ntau = elmN(itaurine,:);
Cgly = elmC(iglycine,:); Ngly = elmN(iglycine,:);

elmC = sum(elmC,1);elmN = sum(elmN,1);
[G, ID] = findgroups(rI.Species);
MeanC = splitapply(@mean,elmC,G'); MeanC = MeanC([4,2,1,3]);
MeanCt = splitapply(@mean,Ctau,G'); MeanCt = MeanCt([4,2,1,3])./1e6;
MeanNt = splitapply(@mean,Ntau,G'); MeanNt = MeanNt([4,2,1,3])./1e6;
MeanCg = splitapply(@mean,Cgly,G'); MeanCg = MeanCg([4,2,1,3])./1e6;
MeanNg = splitapply(@mean,Ngly,G'); MeanNg = MeanNg([4,2,1,3])./1e6;
MeanN = splitapply(@mean,elmN,G'); MeanN = MeanN([4,2,1,3]);
MeanC = MeanC./1e6;MeanN = MeanN./1e6;

CT = DCN{1,2:end} - DCN{1,1};
NT = DCN{2,2:end} - DCN{2,1};
DecoyBar = CT - MeanC;
DecoyBarN = NT - MeanN;

subplot(1,2,1)
b = bar([MeanC',DecoyBar'], "stacked");
b(1).FaceColor = CP1{1};
b(2).FaceColor = CP1{2};
ax = gca;
order = ["Control","Px", "Cp", "Am", "Eu"];
xticklabels(order(2:5))
ylabel("DOC Relative to Control, \muM")
set(ax, "Box", "off")
percents = string(round(100.*MeanC./CT,1))+ "%";
text(b(1).XEndPoints, b(1).YEndPoints, percents, 'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

subplot(1,2,2)
bn = bar([MeanN',DecoyBarN'], "stacked");
bn(1).FaceColor = CP1{1};
bn(2).FaceColor = CP1{2};
ax = gca;
xticklabels(order(2:5))
ylabel("TDN Relative to Control, \muM")
set(ax, "Box", "off")
legend({"Metabolites","Other"})
percentsN = string(round(100.*MeanN./NT,1))+ "%";
text(bn(1).XEndPoints, bn(1).YEndPoints, percentsN, 'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

%% Bar Plot of Excretion rates. 
mtabsofinterest = {'4-aminobenzoic acid','citrulline','taurine','homoserine betaine','isethionate'}; %,'DHPS','arginine','proline','tryptophan','putrescine'};
[~,ia] = ismember(mtabsofinterest,mtabNames);

GraphMtabs = mtabData_pmol(ia,LabelOrder)';
figure
bar(GraphMtabs)
ax = gca;
ax.YScale = 'log';
ax.XTick = 1:22;
legend(mtabsofinterest)

%% Heatmap of Metabolites

HeatMapMtabs = mtabData_pmol(:,LabelOrder);
HeatMapMtabs(HeatMapMtabs==0) = NaN;
%HeatTable = array2table(HeatMapMtabs, "VariableNames",LabelsOrdered + string(1:22)',"RowNames",mtabNames);
h = heatmap(nicenames, LabelsOrdered+ string(1:22)', HeatMapMtabs');
h.ColorLimits = [1,10];
h.ColorScaling = "log";

colors = [CP1{1};CP1{4};CP1{5};CP1{2}];
npts = 1000;
h.Colormap = flip(cmapper(colors,npts));

%% Graphs
if 1
    w = waitbar(0,'','Name','Plotting errorbars...');
    for i = 1:length(mtabNames)
        waitbar(i/length(mtabNames),w, mtabNames(i));
        f1 = figure('Name',mtabNames(i),'Color','none','Position',...
            [500 300 700 500],'Units','inches',"Visible","off");
        ax1 = axes(f1, 'Position', [0.08 0.6 0.41 0.35], 'Units',...
            'normalized');
        ax2 = axes(f1, 'Position', [0.52 0.6 0.41 0.35], 'Units',...
            'normalized');
        ax3 = axes(f1, 'Position', [0.08 0.08 0.41 0.38], 'Units',...
            'normalized');
        ax4 = axes(f1, 'Position', [0.52 0.08 0.41 0.38], 'Units',...
            'normalized');
        
        % The concentration plot.
        
        line1 = plot(ax1,[0,25],[LOD_pM(i),LOD_pM(i)],'LineStyle','--','Color','k');
        hold(ax1)
        line2 = plot(ax1,[0,25],[MaxStd_pM(i),MaxStd_pM(i)],'LineStyle',':','Color','k');
        
        eb1 = errorbar(ax1, meantimes, ctrls(i,:),...
            ste_ctrls(i,:),ste_ctrls(i,:),...
            'Color',CP1{1},'LineWidth',1.5, 'LineStyle', 'none');
        eb2 = errorbar(ax1, Pxtimes, Px_pM(i,:),...
            Px_stde_pM(i,:),Px_stde_pM(i,:),...
            'Color',CP1{3},'LineWidth',1.5, 'LineStyle', 'none');
        eb3 = errorbar(ax1, meantimes(3), Amph_pM(i,:),...
            Amph_stde_pM(i,:),Amph_stde_pM(i,:),...
            'Color',CP1{4},'LineWidth',1.5, 'LineStyle', 'none');
        eb4 = errorbar(ax1, meantimes(3), Clio_pM(i,:),...
            Clio_stde_pM(i,:),Clio_stde_pM(i,:),...
            'Color',CP1{5},'LineWidth',1.5, 'LineStyle', 'none');
        eb5 = errorbar(ax1, meantimes(3), Euph_pM(i,:),...
            Euph_stde_pM(i,:),Euph_stde_pM(i,:),...
            'Color',CP1{2},'LineWidth',1.5, 'LineStyle', 'none');
        ylim([0,max([1,Px_pM(i,:),Amph_pM(i,:),Clio_pM(i,:),Euph_pM(i,:)])])
      
        
        % The legend.
        legend(ax1, {'LowLim','HighLim','ctrl','{\it P. xiphias}','{\it Amphipoda spp.}',...
            '{\it C. pyrimidata}', '{\it Euphasiid spp.}'}, 'Orientation', 'horizontal',...
            'Position', [0.1 0.48 0.82 0.04])
        
        
        % The control subtraction.
        eb6 = errorbar(ax2, meantimes, ctrls(i,:)-meanctrl0(i),...
            ste_ctrls(i,:),ste_ctrls(i,:),...
            'Color',CP1{1},'LineWidth',1.5, 'LineStyle', 'none');
        hold(ax2)
        eb7 = errorbar(ax2, Pxtimes, Px_pM_subctrl(i,:),...
            Px_stde_pM(i,:),Px_stde_pM(i,:),...
            'Color',CP1{3},'LineWidth',1.5, 'LineStyle', 'none');
        eb8 = errorbar(ax2, meantimes(3), Amph_pM_subctrl(i,:),...
            Amph_stde_pM(i,:),Amph_stde_pM(i,:),...
            'Color',CP1{4},'LineWidth',1.5, 'LineStyle', 'none');
        eb9 = errorbar(ax2, meantimes(3), Clio_pM_subctrl(i,:),...
            Clio_stde_pM(i,:),Clio_stde_pM(i,:),...
            'Color',CP1{5},'LineWidth',1.5, 'LineStyle', 'none');
        eb10 = errorbar(ax2, meantimes(3), Euph_pM_subctrl(i,:),...
            Euph_stde_pM(i,:),Euph_stde_pM(i,:),...
            'Color',CP1{2},'LineWidth',1.5, 'LineStyle', 'none'); 
        ylim([0,max([1,Px_pM_subctrl(i,:),Amph_pM_subctrl(i,:),Clio_pM_subctrl(i,:),Euph_pM_subctrl(i,:)])])

        % The inventory.
        eb11 = errorbar(ax3, meantimes, (ctrls(i,:)-meanctrl0(i)).*0.06,...
            ste_ctrls(i,:).*0.06,ste_ctrls(i,:).*0.06,...
            'Color',CP1{1},'LineWidth',1.5, 'LineStyle', 'none');
        hold(ax3)
        eb12 = errorbar(ax3, Pxtimes, Px_pmol(i,:),...
            Px_stde_pmol(i,:),Px_stde_pmol(i,:),...
            'Color',CP1{3},'LineWidth',1.5, 'LineStyle', 'none');
        eb13 = errorbar(ax3, meantimes(3), Amph_pmol(i,:),...
            Amph_stde_pmol(i,:),Amph_stde_pmol(i,:),...
            'Color',CP1{4},'LineWidth',1.5, 'LineStyle', 'none');
        eb14 = errorbar(ax3, meantimes(3), Clio_pmol(i,:),...
            Clio_stde_pmol(i,:),Clio_stde_pmol(i,:),...
            'Color',CP1{5},'LineWidth',1.5, 'LineStyle', 'none');
        eb15 = errorbar(ax3, meantimes(3), Euph_pmol(i,:),...
            Euph_stde_pmol(i,:),Euph_stde_pmol(i,:),...
            'Color',CP1{2},'LineWidth',1.5, 'LineStyle', 'none');
        ylim([0,max([1,Px_pmol(i,:),Amph_pmol(i,:),Clio_pmol(i,:),Euph_pmol(i,:)])])
        
        % The mass-normalized inventory.
        eb16 = errorbar(ax4, meantimes, ctrls(i,:),...
            ste_ctrls(i,:),ste_ctrls(i,:),...
            'Color',CP1{1},'LineWidth',1.5, 'LineStyle', 'none', 'Visible', 'off');
        hold(ax4)
        eb17 = errorbar(ax4, Pxtimes, Px_pmol_mgdry_hr(i,:),...
            Px_stde_pmol_mgdry_hr(i,:),Px_stde_pmol_mgdry_hr(i,:),...
            'Color',CP1{3},'LineWidth',1.5, 'LineStyle', 'none');
        eb18 = errorbar(ax4, meantimes(3), Amph_pmol_mgdry_hr(i,:),...
            Amph_stde_pmol_mgdry_hr(i,:),Amph_stde_pmol_mgdry_hr(i,:),...
            'Color',CP1{4},'LineWidth',1.5, 'LineStyle', 'none');
        eb19 = errorbar(ax4, meantimes(3), Clio_pmol_mgdry_hr(i,:),...
            Clio_stde_pmol_mgdry_hr(i,:),Clio_stde_pmol_mgdry_hr(i,:),...
            'Color',CP1{5},'LineWidth',1.5, 'LineStyle', 'none');
        eb20 = errorbar(ax4, meantimes(3), Euph_pmol_mgdry_hr(i,:),...
            Euph_stde_pmol_mgdry_hr(i,:),Euph_stde_pmol_mgdry_hr(i,:),...
            'Color',CP1{2},'LineWidth',1.5, 'LineStyle', 'none');
        ylim([0,max([1,Px_pmol_mgdry_hr(i,:),Amph_pmol_mgdry_hr(i,:),Clio_pmol_mgdry_hr(i,:),Euph_pmol_mgdry_hr(i,:)])])
        
        ax1.XLabel.String = 'time, h';
        ax1.YLabel.String =  mtabNames(i)+' pM';
        ax2.XLabel.String = 'time, h';
        ax2.YLabel.String = '\DeltaC (-t_o control), pM';
        ax3.XLabel.String = 'time, h';
        ax3.YLabel.String = '\DeltaC*V (inventory), pmol';
        ax4.XLabel.String = 'time, h';
        ax4.YLabel.String = 'pmol mg^{-1} h^{-1}';
        
        set(ax1, 'XLim', [0, 13])
        set(ax2, 'XLim', [0, 13],'YAxisLocation', 'right')
        set(ax3, 'XLim', [0, 13])
        set(ax4, 'XLim', [0, 13],'YAxisLocation', 'right')
        
        if 0
            break
        end
        
        name = string(outdir) + '/' + "PeeRatesMatchedCtrl" + '.pdf';
        exportgraphics(f1, name, 'Append',  true)
        close all
        clear f1 ax1 ax2 ax3 ax4 eb1 eb2 eb3 eb4 eb5 eb6 eb7 eb8 eb9 eb10 eb11 eb12 eb13 eb14 eb15 eb16 line1 line2
    end
end

%% Summary Figure. 

%goodStuff = ratingFlags<3;
reducedRatePx = Px_pmol_mgdry_hr;%(goodStuff);
reducedErrPx = Px_stde_pmol_mgdry_hr;%(goodStuff);
reducedRateAmph = Amph_pmol_mgdry_hr;%(goodStuff);
reducedRateEuph = Euph_pmol_mgdry_hr;
reducedErrAmph = Amph_stde_pmol_mgdry_hr;%(goodStuff);
reducedErrEuph = Euph_stde_pmol_mgdry_hr;
reducedRateClio = Clio_pmol_mgdry_hr;%(goodStuff);
reducedErrClio = Clio_stde_pmol_mgdry_hr;%(goodStuff);

reducedNames = scattertabs;%(goodStuff);
elem = mtabElem; 
%elem = elem(goodStuff,:);

% Trim anything <0
reducedRates = [reducedRatePx(:,2), reducedRateAmph, reducedRateClio, reducedRateEuph];
reducedRates(reducedRates<0) = 0;
reducedErr = [reducedErrPx(:,2), reducedErrAmph, reducedErrClio, reducedRateEuph];
ibad = nansum(reducedRates,2) == 0 |...
    isnan(nansum(reducedRates,2));
reducedRates(ibad,:) = [];
reducedNames(ibad,:) = [];
reducedErr(ibad,:)= [];
elem(ibad,:)=[];


if 1
    f1 = figure('Name',"AllEstimates",'Color','none','Position',...
        [50 30 800 1000],'Units','inches');
    
    b1 = bar(categorical(reducedNames(max(reducedRates>1,[],2))),...
        reducedRates(max(reducedRates>1,[],2),:), 'FaceColor','flat');
    % b1(1).XData = categorical(reducedNames(max(reducedRates>1,[],2)));
    % b1(2).XData = categorical(reducedNames(max(reducedRates>1,[],2)));
    % b1(3).XData = categorical(reducedNames(max(reducedRates>1,[],2)));
    % b1(4).XData = categorical(reducedNames(max(reducedRates>1,[],2)));
    b1(1).CData = CP1{2};
    b1(2).CData = repelem(CP1{4},length(reducedNames(max(reducedRates>1,[],2))),1);
    b1(3).CData = repelem(CP1{5},length(reducedNames(max(reducedRates>1,[],2))),1);
    b1(4).CData = repelem(CP1{3},length(reducedNames(max(reducedRates>1,[],2))),1);
    set(gca, 'XTickLabel', categorical(reducedNames(max(reducedRates>1,[],2))))
    set(gca,... 'XTick', 1:length(reducedNames(max(reducedRates>1,[],2))),...
        'XTickLabelRotation', 90, 'YScale', 'log', 'TickLength', [0 0], ...
        'YGrid', 'on', 'YLim', [1, max(max(reducedRates(reducedRates>1)))+10]);
    ylabel('pmol mg^{-1} h^{-1}', 'Interpreter', 'tex')
    legend(gca, {'{\it P. xiphias}','{\it Amphipoda spp.}',...
        '{\it C. cuspidata}', '{\it Euphasiid spp.}'}, 'Orientation', 'vertical','Location',...
        'northwest', 'FontSize', 8)
    saveas(f1, outdir + "/"+"EstimatedRates.pdf")
end

%% Sorting Rate Data: Eliminate Above-Curve and High-Error.

i12h = (sI.Nominal_Duration_h==12);
iDead = sI.Notes=="DEAD";

% Flag 1 is 1 if the sample is above the standard curve. 
Pxf1 = mtabData_pM_Reorder(:,i12h & iPx&~iDead)>MaxStd_pM;
Cliof1 = mtabData_pM_Reorder(:,i12h & iClio&~iDead)>MaxStd_pM;
Amphf1 = mtabData_pM_Reorder(:,i12h & iAmph&~iDead)>MaxStd_pM;
Euphf1 = mtabData_pM_Reorder(:,i12h & iEuph&~iDead)>MaxStd_pM;

maxoutTable = table(mtabNames, mean(Pxf1,2), mean(Cliof1,2),...
    mean(Amphf1,2), mean(Euphf1,2),...
    'VariableNames',{'Mtab','Px','Clio','Amph','Euph'});

% Flag 2 is 1 if the sample is below the LOD.
Pxf2 = mtabData_pM_Reorder(:,i12h & iPx&~iDead)<LOD_pM;
Cliof2 = mtabData_pM_Reorder(:,i12h & iClio&~iDead)<LOD_pM;
Amphf2 = mtabData_pM_Reorder(:,i12h & iAmph&~iDead)<LOD_pM;
Euphf2 = mtabData_pM_Reorder(:,i12h & iEuph&~iDead)<LOD_pM;

LODTable = table(mtabNames, mean(Pxf2,2), mean(Cliof2,2),...
    mean(Amphf2,2), mean(Euphf2,2),...
    'VariableNames',{'Mtab','Px','Clio','Amph','Euph'});

% Flag 3 is 1 if the samples at t12 are not significantly different from
% zero.
Pxf3 = ttest(mtabData_pM_Reorder(:,i12h & iPx&~iDead)')'; Pxf3(isnan(Pxf3))=0;
Pxf3 = ~Pxf3;
Cliof3 = ones(size(mtabData_pM_Reorder(:,i12h & iClio&~iDead)));%ONLY ONE LIVE SAMPLE ttest(mtabData_pM_Reorder(:,i12h & iClio&~iDead)')'; 
Cliof3(isnan(Cliof3))=0;
Cliof3 = ~Cliof3;
Amphf3 = ttest(mtabData_pM_Reorder(:,i12h & iAmph&~iDead)')'; Amphf3(isnan(Amphf3))=0;
Amphf3 = ~Amphf3;
Euphf3 = ttest(mtabData_pM_Reorder(:,i12h & iEuph&~iDead)')'; Euphf3(isnan(Euphf3))=0;
Euphf3 = ~Euphf3;

ttestTable = table(mtabNames, mean(Pxf3,2), mean(Cliof3,2),...
    mean(Amphf3,2), mean(Euphf3,2),...
    'VariableNames',{'Mtab','Px','Clio','Amph','Euph'});

% Flag 4 is 1 if the samples at t12 are not significantly different from
% the time-matched controls. 
Pxf4 = ttest2(mtabData_pM_Reorder(:,i12h & iPx&~iDead)',t12ctrl')'; Pxf4(isnan(Pxf4))=0;
Pxf4 = ~Pxf4;
Cliof4 = ttest2(t12ctrl', mtabData_pM_Reorder(:,i12h & iClio&~iDead)')'; Cliof4(isnan(Cliof4))=0;
Cliof4 = ~Cliof4;
Amphf4 = ttest2(mtabData_pM_Reorder(:,i12h & iAmph&~iDead)',t12ctrl')'; Amphf4(isnan(Amphf4))=0;
Amphf4 = ~Amphf4;
Euphf4 = ttest2(mtabData_pM_Reorder(:,i12h & iEuph&~iDead)',t12ctrl')'; Euphf4(isnan(Euphf4))=0;
Euphf4 = ~Euphf4;

ttestctrlTable = table(mtabNames, mean(Pxf4,2), mean(Cliof4,2),...
    mean(Amphf4,2), mean(Euphf4,2),...
    'VariableNames',{'Mtab','Px','Clio','Amph','Euph'});

%% What are the main disqualifiers? Do I need to dilute and rerun?

AboveMax = table(mtabNames);

AboveMax.FracPx = sum(Pxf1, 2)./2;
AboveMax.FracClio = sum(Cliof1,2)./2;
AboveMax.FracAmph = sum(Amphf1,2)./3;
AboveMax.FracEuph = sum(Euphf1,2)./3;
AboveMax.Avg = sum([Pxf1,Cliof1,Amphf1,Euphf1],2)./10;
AboveMax.All = sum([Pxf1,Cliof1,Amphf1,Euphf1],2);

writetable(AboveMax, "../datasets/AboveMax.csv")

AboveMax_noProb = table(mtabNames);

OnlyAbovePx = Pxf1 & ~(Pxf2 | Pxf3 | Pxf4);
OnlyAboveClio = Cliof1 & ~(Cliof2 | Cliof3 | Cliof4);
OnlyAboveEuph = Euphf1 & ~(Euphf2 | Euphf3 | Euphf4);
OnlyAboveAmph = Amphf1 & ~(Amphf2 | Amphf3 | Amphf4);

AboveMax_noProb.FracPx = sum(OnlyAbovePx, 2)./2;
AboveMax_noProb.FracClio = sum(OnlyAboveClio, 2)./2;
AboveMax_noProb.FracAmph = sum(OnlyAboveAmph, 2)./3;
AboveMax_noProb.FracEuph = sum(OnlyAboveEuph, 2)./3;
AboveMax_noProb.Avg = sum([OnlyAbovePx,OnlyAboveClio,OnlyAboveAmph,OnlyAboveEuph],2)./10;
AboveMax_noProb.All = sum([OnlyAbovePx,OnlyAboveClio,OnlyAboveAmph,OnlyAboveEuph],2);

writetable(AboveMax_noProb, "../datasets/AboveMax_noProb.csv")

%% Creating hyper-reduced versions of the rate data using the flags. 

iReject_Px = Pxf1 | Pxf2 | Pxf3 | Pxf4;
iReject_Clio = Cliof1 | Cliof2 | Cliof3 | Cliof4;
iReject_Amph = Amphf1 | Amphf2 | Amphf3 | Amphf4;
iReject_Euph = Euphf1 | Euphf2 | Euphf3 | Euphf4;

% The inventory
Px_pmol0 = mtabData_pmol(:,t12i&iPx&~iDead); Px_pmol0(iReject_Px) = NaN;
Px_pmol = mean(Px_pmol0,2,'omitnan');
Px_stde_pmol = std(Px_pmol0,[],2,'omitnan');
clear Px_pmol0

Amph_pmol0 = mtabData_pmol(:,t12i&iAmph&~iDead); Amph_pmol0(iReject_Amph) = NaN;
Amph_pmol = mean(Amph_pmol0,2,'omitnan');
Amph_stde_pmol = std(Amph_pmol0,[],2,'omitnan');
clear Amph_pmol0

Clio_pmol0 = mtabData_pmol(:,t12i&iClio&~iDead); Clio_pmol0(iReject_Clio) = NaN;
Clio_pmol = mean(Clio_pmol0,2,'omitnan');
Clio_stde_pmol = std(Clio_pmol0,[],2,'omitnan');
clear Clio_pmol0

Euph_pmol0 = mtabData_pmol(:,t12i&iEuph&~iDead); Euph_pmol0(iReject_Euph) = NaN;
Euph_pmol = mean(Euph_pmol0,2,'omitnan');
Euph_stde_pmol = std(Euph_pmol0,[],2,'omitnan');
clear Euph_pmol0

% Now, the time-normalized inventory.
Px_pmol_mgdry_hr0 = mtabData_pmol_mgdry_hr(:,t12i&iPx&~iDead); Px_pmol_mgdry_hr0(iReject_Px) = NaN;
Px_pmol_mgdry_hr = mean(Px_pmol_mgdry_hr0,2,'omitnan');
Px_stde_pmol_mgdry_hr = std(Px_pmol_mgdry_hr0,[],2,'omitnan');
clear Px_pmol_mgdry_hr0

Amph_pmol_mgdry_hr0 = mtabData_pmol_mgdry_hr(:,t12i&iAmph&~iDead); Amph_pmol_mgdry_hr0(iReject_Amph) = NaN;
Amph_pmol_mgdry_hr = mean(Amph_pmol_mgdry_hr0,2,'omitnan');
Amph_stde_pmol_mgdry_hr = std(Amph_pmol_mgdry_hr0,[],2,'omitnan');
clear Amph_pmol_mgdry_hr0

Clio_pmol_mgdry_hr0 = mtabData_pmol_mgdry_hr(:,t12i&iClio&~iDead); Clio_pmol_mgdry_hr0(iReject_Clio) = NaN;
Clio_pmol_mgdry_hr = mean(Clio_pmol_mgdry_hr0,2,'omitnan');
Clio_stde_pmol_mgdry_hr = std(Clio_pmol_mgdry_hr0,[],2,'omitnan');
clear Clio_pmol_mgdry_hr0

Euph_pmol_mgdry_hr0 = mtabData_pmol_mgdry_hr(:,t12i&iEuph&~iDead); Euph_pmol_mgdry_hr0(iReject_Euph) = NaN;
Euph_pmol_mgdry_hr = mean(Euph_pmol_mgdry_hr0,2,'omitnan');
Euph_stde_pmol_mgdry_hr = std(Euph_pmol_mgdry_hr0,[],2,'omitnan');
clear Euph_pmol_mgdry_hr0


%% Summary Figure Redux

reducedRatePx = Px_pmol_mgdry_hr;
reducedErrPx = Px_stde_pmol_mgdry_hr;
reducedRateAmph = Amph_pmol_mgdry_hr;
reducedRateEuph = Euph_pmol_mgdry_hr;
reducedErrAmph = Amph_stde_pmol_mgdry_hr;
reducedErrEuph = Euph_stde_pmol_mgdry_hr;
reducedRateClio = Clio_pmol_mgdry_hr;
reducedErrClio = Clio_stde_pmol_mgdry_hr;

reducedNames = scattertabs;
elem = mtabElem; 

% Trim anything <0
reducedRates = [reducedRatePx, reducedRateAmph, reducedRateClio, reducedRateEuph];
reducedRates(reducedRates<0) = 0;
reducedErr = [reducedErrPx, reducedErrAmph, reducedErrClio, reducedRateEuph];
ibad = nansum(reducedRates,2) == 0 |...
    isnan(nansum(reducedRates,2));
reducedRates(ibad,:) = [];
reducedNames(ibad,:) = [];
reducedErr(ibad,:)= [];
elem(ibad,:)=[];


if 1
    f1 = figure('Name',"AllEstimates_trimmed",'Color','none','Position',...
        [50 30 800 1000],'Units','inches');
    
    b1 = bar(categorical(reducedNames(max(reducedRates>1,[],2))),...
        reducedRates(max(reducedRates>1,[],2),:), 'FaceColor','flat');
    b1(1).CData = CP1{2};
    b1(2).CData = repelem(CP1{4},length(reducedNames(max(reducedRates>1,[],2))),1);
    b1(3).CData = repelem(CP1{5},length(reducedNames(max(reducedRates>1,[],2))),1);
    b1(4).CData = repelem(CP1{3},length(reducedNames(max(reducedRates>1,[],2))),1);
    set(gca, 'XTickLabel', categorical(reducedNames(max(reducedRates>1,[],2))))
    set(gca,... 'XTick', 1:length(reducedNames(max(reducedRates>1,[],2))),...
        'XTickLabelRotation', 90, 'YScale', 'log', 'TickLength', [0 0], ...
        'YGrid', 'on', 'YLim', [1, max(max(reducedRates(reducedRates>1)))+10]);
    ylabel('pmol mg^{-1} h^{-1}', 'Interpreter', 'tex')
    legend(gca, {'{\it P. xiphias}','{\it Amphipoda spp.}',...
        '{\it C. cuspidata}', '{\it Euphasiid spp.}'}, 'Orientation', 'vertical','Location',...
        'northwest', 'FontSize', 8)
    saveas(f1, outdir + "/"+"EstimatedRates_trimmed.pdf")
end
%% Analysis of Rate Data. 

rateEst = table();
rateEst.mtabName = mtabNames;

rateEst.Px = nanmean(mtabData_pmol_mgdry_hr(:,i12h & iPx),2);
rateEst.Clio = nanmean(mtabData_pmol_mgdry_hr(:,i12h & iClio),2);
rateEst.Amph = nanmean(mtabData_pmol_mgdry_hr(:,i12h & iAmph),2);
rateEst.Euph = nanmean(mtabData_pmol_mgdry_hr(:,i12h & iEuph),2);
rateEst.PxFlag = nanmean(mtabData_pM_Reorder(:,i12h & iPx),2)>MaxStd_pM;
rateEst.ClioFlag = nanmean(mtabData_pM_Reorder(:,i12h & iClio),2)>MaxStd_pM;
rateEst.AmphFlag = nanmean(mtabData_pM_Reorder(:,i12h & iAmph),2)>MaxStd_pM;
rateEst.EuphFlag = nanmean(mtabData_pM_Reorder(:,i12h & iEuph),2)>MaxStd_pM;
rateEst.ErrPx = Px_stde_pmol_mgdry_hr;% (:,2);
rateEst.ErrAmph = Amph_stde_pmol_mgdry_hr;
rateEst.ErrClio = Clio_stde_pmol_mgdry_hr;
rateEst.ErrEuph = Euph_stde_pmol_mgdry_hr;
rateEst.Px(rateEst.Px<=0 | isnan(rateEst.Px))=0;
rateEst.Amph(rateEst.Amph<=0 | isnan(rateEst.Amph))=0;
rateEst.Clio(rateEst.Clio<=0 | isnan(rateEst.Clio))=0;
rateEst.Euph(rateEst.Euph<=0 | isnan(rateEst.Euph))=0;

rateEst.General = nanmean(mtabData_pmol_mgdry_hr(:,i12h & (iClio| iAmph | iPx |iEuph)),2);
rateEst.C_rate_Px = sum([rateEst.Px],2).*mtabElem.C;
rateEst.N_rate_Px = sum([rateEst.Px],2).*mtabElem.N;
rateEst.O_rate_Px = sum([rateEst.Px],2).*mtabElem.O;
rateEst.C_err_Px = sum([rateEst.ErrPx],2).*mtabElem.C;
rateEst.N_err_Px = sum([rateEst.ErrPx],2).*mtabElem.N;
rateEst.O_err_Px = sum([rateEst.ErrPx],2).*mtabElem.O;
rateEst.C_rate_Amph = sum([rateEst.Amph],2).*mtabElem.C;
rateEst.N_rate_Amph = sum([rateEst.Amph],2).*mtabElem.N;
rateEst.O_rate_Amph = sum([rateEst.Amph],2).*mtabElem.O;
rateEst.C_err_Amph = sum([rateEst.ErrAmph],2).*mtabElem.C;
rateEst.N_err_Amph = sum([rateEst.ErrAmph],2).*mtabElem.N;
rateEst.O_err_Amph = sum([rateEst.ErrAmph],2).*mtabElem.O;
rateEst.C_rate_Clio = sum([rateEst.Clio],2).*mtabElem.C;
rateEst.N_rate_Clio = sum([rateEst.Clio],2).*mtabElem.N;
rateEst.O_rate_Clio = sum([rateEst.Clio],2).*mtabElem.O;
rateEst.C_err_Clio = sum([rateEst.ErrClio],2).*mtabElem.C;
rateEst.N_err_Clio = sum([rateEst.ErrClio],2).*mtabElem.N;
rateEst.O_err_Clio = sum([rateEst.ErrClio],2).*mtabElem.O;
rateEst.C_rate_Euph = sum([rateEst.Euph],2).*mtabElem.C;
rateEst.N_rate_Euph = sum([rateEst.Euph],2).*mtabElem.N;
rateEst.O_rate_Euph = sum([rateEst.Euph],2).*mtabElem.O;
rateEst.C_err_Euph = sum([rateEst.ErrEuph],2).*mtabElem.C;
rateEst.N_err_Euph = sum([rateEst.ErrEuph],2).*mtabElem.N;
rateEst.O_err_Euph = sum([rateEst.ErrEuph],2).*mtabElem.O;
rateEst(rateEst.General <=0 | isnan(rateEst.General),:)=[];
elem(rateEst.General <=0 | isnan(rateEst.General),:)=[];


writetable(rateEst, '../datasets/rateEstimates_22Sept.xlsx')
