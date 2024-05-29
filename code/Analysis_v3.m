% Noah Germolus 06 Feb 2024
% This is v3 of Analysis. It will remove dead animals from processing, as
% well as calculate the rates and regressions differently.

clear; clc;  
close all

% Loads data and metadata.
loadstuff

outdir = '../figs/zoopfigs_v3';
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

% Information for constructing a color gradient later on. 
colors = [CP1{2};CP1{1};CP1{3}];
col2 = [CP1{4}; CP1{3}];
npts = 1000;

%% Calculations.
%
%   I want to calculate several things. 
%   -change in concentration: THIS IS DIFFERENT IN v3. Now, the difference
%   between the time-matched control and the t0 control is subtracted from
%   the time-point. This may result in higher numbers for many samples. 
%   -change in inventory: same calculation, using the volume-normalized
%       mtabData_pmol
%   -change at 12h on a dry-mass basis. This will utilize the
%       volume-normalized inventories, but additionally divided by the dry
%       weight of the bugs. 
%

% Dead animals
idead = sI.Notes == "DEAD";
sI(idead,:) = [];
% Determine the actual incubation durations. 
sI.duration = sI.TimeStop - sI.TimeStart;

% Create useful indices for controls.
t0ctrl_i = find(sI.Nominal_Duration_h == 0);
t6ctrl_i = find(sI.Nominal_Duration_h == 6 & sI.Species == "CTRL");
t12ctrl_i = find(sI.Nominal_Duration_h == 12 & sI.Species == "CTRL");
iallctrl = sI.Species=="CTRL"; 

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
% Get rid of deads
mtabData_pM_Reorder = mtabData_pM_Reorder(:,~idead);
% Subrtract time-matched controls.
mtabData_pM_Reorder_subctrl = mtabData_pM_Reorder;
mtabData_pM_Reorder_subctrl(:,sI.Nominal_Duration_h==0) = mtabData_pM_Reorder_subctrl(:,sI.Nominal_Duration_h==0) - meanctrl0;
mtabData_pM_Reorder_subctrl(:,sI.Nominal_Duration_h==6) = mtabData_pM_Reorder_subctrl(:,sI.Nominal_Duration_h==6) - (meanctrl6-meanctrl0);
mtabData_pM_Reorder_subctrl(:,sI.Nominal_Duration_h==12) = mtabData_pM_Reorder_subctrl(:,sI.Nominal_Duration_h==12) - (meanctrl12-meanctrl0);
mtabData_pM_Reorder_subctrl(mtabData_pM_Reorder_subctrl<0) = NaN;

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

%%  Creating a screening process for metabolites that are too high in the
% control or else just bad ones that slipped past the automated QC

mtab_pM_ctrl = mtabData_pM_Reorder(:,iallctrl);
mtab_pM_notctrl = mtabData_pM_Reorder(:,~iallctrl);
mremove = (10.*mean(mtab_pM_ctrl,2,"omitnan")>mean(mtab_pM_notctrl,2,"omitnan"));
rmNames = ["2'deoxyguanosine Na","S-(5'-adenosyl)-L-homocysteine",...
     "adenosine", "serine 2", "cysteine dimer",...
     "ornithine 2"]'; %"4-aminobenzoic acid",
[~, ibad] = ismember(rmNames, mtabNames); mremove(ibad) = 1;

clear rmNames mtab_pM_notctrl mtab_pM_ctrl

%% Just giving you some idea of what the top ten
% mean excreted compounds are for P. xiphias, 6 h
plot_top10

%% Heatmap plot
% Originally, this section contained a lot about ordinating and clustering
% the data based on Bray-Curtis dissimilarity. This was not a particularly
% useful metric, and while that analysis still exists in Analysis_OneMode,
% I will not reproduce it here, instead favoring an approach that simply
% groups the species based on taxonomy. 

% I created the variable TaxOrder which will order the species by
% high-order (phylum to order) taxonomy: CTRL, pteropod, copepod, amphipod,
% euphausiids
[~,LabelOrder] = sort(sI.TaxOrder(~iallctrl));
LabelsOrdered = sI.AccurateSpecies(~iallctrl);
LabelsOrdered = LabelsOrdered(LabelOrder);
HeatMapMtabs = mtabData_pmol_mgdry_hr(~mremove, ~iallctrl);
HeatMapMtabs = round(HeatMapMtabs(:,LabelOrder));
HeatMapMtabs(HeatMapMtabs<=0) = NaN;
ihrm = (sum(isnan(HeatMapMtabs),2)==size(HeatMapMtabs,2));
% LiveDead = sI.Notes(~iallctrl); LiveDead = string(LiveDead(LabelOrder));
TimePointInfo = sI.Nominal_Duration_h(~iallctrl); 
TimePointInfo = string(TimePointInfo(LabelOrder));
hxdata = nicenames(~mremove);
hydata = LabelsOrdered+ " " + TimePointInfo +" h "+...% + LiveDead+ " " +...
    string(1:size(TimePointInfo,1))';
hxdata(ihrm,:)=[];
HeatMapMtabs(ihrm,:)=[];
% Make the heatmap.
h = heatmap(hydata, hxdata, HeatMapMtabs);
% h.ColorLimits = [1,10];
h.ColorScaling = "scaledcolumns";
h.Colormap = flip(cmapper(col2,npts));
h.FontSize = 12;
h.ColorbarVisible = "off";
axp = struct(gca);       %you will get a warning
axp.Axes.XAxisLocation = 'top';

% Add the "Above the standard curve" markers
if 0
    MaxRM = MaxStd_pM(~mremove);
    AboveMax = mtabData_pM_Reorder(~mremove,~iallctrl);
    AboveMax = AboveMax(:,LabelOrder)>MaxRM;
    AboveMax = flip(AboveMax',1);
    ax = axes(gcf, "Position", h.Position, "Units","normalized","Color", "none");
    ax.Box = "off";
    set(ax, "XTick", [], "YTick", [])
    ax.XLim = [0 size(AboveMax,2)]; ax.YLim = [0 size(AboveMax,1)];
    [X,Y] = meshgrid(0:size(AboveMax,2),0:size(AboveMax,1));
    t = text(X(AboveMax), Y(AboveMax),"*");
end

%% Heatmap plot: reduced to ten.
% Originally, this section contained a lot about ordinating and clustering
% the data based on Bray-Curtis dissimilarity. This was not a particularly
% useful metric, and while that analysis still exists in Analysis_OneMode,
% I will not reproduce it here, instead favoring an approach that simply
% groups the species based on taxonomy. 

% I created the variable TaxOrder which will order the species by
% high-order (phylum to order) taxonomy: CTRL, pteropod, copepod, amphipod,
% euphausiids

% This new version only uses the top ten overall excreted mtabs. 
[~,LabelOrder] = sort(sI.TaxOrder(~iallctrl));
LabelsOrdered = sI.AccurateSpecies(~iallctrl);
LabelsOrdered = LabelsOrdered(LabelOrder);
HeatMapMtabs = mtabData_pmol_mgdry_hr(~mremove, ~iallctrl);
HeatMapMtabs = round(HeatMapMtabs(:,LabelOrder));
HeatMapMtabs(HeatMapMtabs<=0) = NaN;
ihrm = (sum(isnan(HeatMapMtabs),2)==size(HeatMapMtabs,2));
%LiveDead = sI.Notes(~iallctrl); LiveDead = string(LiveDead(LabelOrder));
TimePointInfo = sI.Nominal_Duration_h(~iallctrl); 
TimePointInfo = string(TimePointInfo(LabelOrder));
hxdata = nicenames(~mremove);
hydata = LabelsOrdered+ " " + TimePointInfo +" h "+... % + LiveDead+ " " +...
    string(1:size(TimePointInfo,1))';
hxdata(ihrm,:)=[];
HeatMapMtabs(ihrm,:)=[];

HeatMapMeds = HeatMapMtabs;
HeatMapMeds(isnan(HeatMapMeds)) = 0;
medians = median(HeatMapMeds,2); 
[~,iam] = sort(medians, "descend");
HeatMapMtabs = HeatMapMtabs(iam(1:10),:);
hxdata = hxdata(iam(1:10));
tooHigh = mtabData_pM_Reorder>MaxStd_pM; 
tooHigh_map = tooHigh(~mremove, ~iallctrl);
tooHigh_map = tooHigh_map(:,LabelOrder);
tooHigh_map(ihrm,:) = [];
tooHigh_map = tooHigh_map(iam(1:10),:);

% Make the heatmap.
f = figure;
h = heatmap(hydata, hxdata, HeatMapMtabs);
% h.ColorLimits = [1,10];
h.ColorScaling = "scaledcolumns";
h.Colormap = flip(cmapper(col2,npts));
h.FontSize = 12;
h.ColorbarVisible = "off";
axp = struct(gca);       %you will get a warning
axp.Axes.XAxisLocation = 'top';
ax = gca;

% Add the "Above the standard curve" markers
if 0
    MaxRM = MaxStd_pM(~mremove);
    AboveMax = mtabData_pM_Reorder(~mremove,~iallctrl);
    AboveMax = AboveMax(:,LabelOrder)>MaxRM;
    AboveMax = flip(AboveMax',1);
    ax2 = axes(f,"Position", ax.Position,"Color", "none");
    ax2.Box = "off";
    set(ax2, "XTick", [], "YTick", [])
    ax2.XLim = [0 size(AboveMax,2)]; ax2.YLim = [0 size(AboveMax,1)];
    [X,Y] = meshgrid(0:size(AboveMax,2),0:size(AboveMax,1));
    t = text(X(AboveMax), Y(AboveMax),"*");
end


%% Bray-Curtis Dissimilarity and ANOSIM
% I am going to conduct this analysis using the codes from the Fathom
% toolbox: f_braycurtis, f_anosim, and f_anosim2
% also requires dependencies, so really you should download the fathom
% toolbox and include it as below.
addpath("F:/Noah Germolus/Documents/MATLAB/Fathom/")

% First-things first is constructing two dissimilarity matrices:
% One that includes all samples on a mole inventory basis...
mn = mtabData_pM_Reorder; mn(isnan(mn))=0; % set NaNs to zero
diss_inv = f_braycurtis(mn); % Generate a distance matrix

% ...and one based on the normalized rates, which eliminates the controls
% through the earlier subtraction step.
mn = mtabData_pmol_mgdry_hr(:,~iallctrl); mn(isnan(mn))=0; % set NaNs to zero
diss_rate = f_braycurtis(mn); % Generate a distance matrix

clear mn

% Next, we do two different ANOSIMs:
% The first, using the diss_inv to see if the controls are significantly
% different from the samples...
grps = iallctrl+1;
%disp(grpn_inv)
ANOSIM_inv = f_anosim(diss_inv,grps,1,1000,1,1);

% clear grps
% % ...and the second, to drill down into the rates themselves without the
% % controls. 
% [grps, grpn_rate] = findgroups(sI.AccurateSpecies(~iallctrl));
% disp(grpn_rate)
% ANOSIM_rate = f_anosim(diss_rate,grps,1,1000,1,1);

clear grps

%% 27 Dec 2023: Different scaling relationship?
% Correlations between mass and 12 h export
in = (sI.Species ~="CTRL");

% Testing here the power-law relationship fit between biomass and
% excretion.
rnonlog = mtabData_pmol(~mremove,in)'./hours(sI.duration(in));
pwrE = log(rnonlog);
pwrE(isinf(pwrE))=NaN;
pwrM = log(sI.dryWeight(in));
% Eliminate things with <3 data points
pwrBad = sum(~(isnan(pwrE)))<3;
pwrE = pwrE(:,~pwrBad); 
pwrNames = nicenames(~mremove);
pwrNames = nicenames(~pwrBad);
rnonlog = rnonlog(:,~pwrBad);
R2 = zeros(size(pwrNames));
p = zeros(size(pwrNames));
m = zeros(size(pwrNames));
b = zeros(size(pwrNames));
n = zeros(size(pwrNames));
resTable = table(pwrNames, m, b, R2, p);
for ii = 1:length(pwrNames)
    lm = fitlm(pwrM,pwrE(:,ii));
    resTable.b(ii) = lm.Coefficients.Estimate(1);
    resTable.m(ii) = lm.Coefficients.Estimate(2);
    resTable.R2(ii) = lm.Rsquared.Adjusted;
    resTable.p(ii) = lm.ModelFitVsNullModel.Pvalue;
    resTable.n(ii) = sum(~isnan(pwrE(:,ii)));
    if lm.Rsquared.Ordinary>0.5
        figure
        x = linspace(min(pwrM),max(pwrM),100);
        y = x.*resTable.m(ii)+resTable.b(ii);   
        plot(exp(x),exp(y))
        hold on 
        scatter(sI.dryWeight(in),rnonlog(:,ii));   
        title([pwrNames(ii)+" R^2= "+string(resTable.R2(ii))]);
        xlabel("Dry Biomass, mg")
        ylabel("R, pmol hr^{-1}", "Interpreter","tex")
        ylim([0,max(pwrE(:,ii))+5])
        ax1 = gca;
        ax1.XScale = "log";
        ax1.YScale = "log";
        % pos = ax1.Position;
        % ax2 = axes("Position",[0.200 0.5 0.3 0.3]);
        % plot(x,y)
        % hold on
        % %plot(x,exp(lm.)
        % scatter(sI.dryWeight(t12an),rnonlog(:,ii));
        % ylim([0,max(rnonlog(rnonlog(:,ii)<max(rnonlog(:,ii)),ii))+5]);
        % xlim([0 4.5])
    end

end
resTable = sortrows(resTable,"R2","descend");

%% Scatterplot; by species

scattertabs = mtabData_pmol_mgdry_hr(~mremove,(t12i|t6i))';
scattertabs(isnan(scattertabs))= 0;
rI = sI(t12i|t6i,:);
ictrl = rI.Species=="CTRL";
rI = rI(~ictrl,:);
scattertabs = scattertabs(~ictrl,:);
scattertabs(scattertabs<0)=0;
% ia12 = rI.Nominal_Duration_h==12&(rI.Species=="PX"|rI.Species=="C. pyrimidata"|rI.Species=="Euph"|rI.Species=="Amphipod (long skinny)");
% clear ia12
kfun = @(m) 10./(6*m);
hline = @(y, x) plot(x, [y,y], "--k", "HandleVisibility","off");

subplot(2,2,1)
% P. xiphias plot.
pxtabs6 = scattertabs(rI.Species == "PX"&rI.Nominal_Duration_h==6,:);
pxtabs12 = scattertabs(rI.Species == "PX"&rI.Nominal_Duration_h==12,:);
meanmass = mean(rI.dryWeight(rI.Species == "PX"));
kthresh = kfun(meanmass);
% below = sum(pxtabs>2.1667,1)==0;
med = median([pxtabs6;pxtabs12],1,"omitnan");
[~,ia] = sort(med, "ascend", "MissingPlacement","last");
% pxtabs6 = pxtabs6(:,ia(~below));
% pxtabs12 = pxtabs12(:,ia(~below));

pxnames = nicenames(~mremove); 
pxnames = pxnames(ia); %pxnames(ia(~below));
pxtabs6 = pxtabs6(:,ia);pxtabs12=pxtabs12(:,ia);
noplot = mean([pxtabs6;pxtabs12],1)==0;
pxnames = pxnames(~noplot);
pxtabs6 = pxtabs6(:,~noplot); pxtabs12 = pxtabs12(:,~noplot);
xPx6 = repmat(1:size(pxnames,1),1,sum(rI.Species=="PX"&rI.Nominal_Duration_h==6))'; %1:sum(~below)
xPx12 = repmat(1:size(pxnames,1),1,sum(rI.Species=="PX"&rI.Nominal_Duration_h==12))'; %1:sum(~below)
bmmPx6 = [pxtabs6(1,:),pxtabs6(2,:),pxtabs6(3,:)]';
bmmPx12 = [pxtabs12(1,:),pxtabs12(2,:),pxtabs12(3,:)]';
sc1 = scatter(xPx6,bmmPx6, 65, CP1{5},"v","LineWidth",2);
hold on
sc2 = scatter(xPx12,bmmPx12, 65, CP1{5},'filled',"^");
ax = gca;
% ax.Color = "none";
% ax.XAxis.Color = soft{1};
% ax.YAxis.Color = soft{1};
% sc1.MarkerEdgeColor = soft{2};
% sc2.MarkerEdgeColor = soft{2};
% sc2.MarkerFaceColor = soft{3};
ax.YLabel.String = "pmol h^{-1} mg^{-1}";
set(ax,"XLim",[0 1+length(pxnames)], "YScale","log");
set(ax,"XTick",1:length(pxnames), "XTickLabels", pxnames,"XTickLabelRotation",45)
h = hline(kthresh, ax.XLim);
%h.Color = soft{2};
%%% FORMAT THESE
% xl6 = log10(bmmPx6); xl6(isinf(xl6))=NaN;
% xl12 = log10(bmmPx12); xl12(isinf(xl12))=NaN;
% lm6 = fitlm(xPx6,xl6); lm12 = fitlm(xPx12,xl12);
% fl6 = plot(xPx6,10.^(predict(lm6,xPx6)), "HandleVisibility","on", "LineStyle","-","Color",[0.5 0.5 0.5], "LineWidth",2.5);
% fl12 = plot(xPx12,10.^(predict(lm12,xPx12)), "HandleVisibility","on", "LineStyle",'-',"Color",[0.2 0.2 0.2], "LineWidth",2.5);
%%%
xL = 1.5:1:(1.5+size(pxnames,1)-2);
xline(xL, "LineStyle","--", "LineWidth",0.5, "Color",.5.*CP1{5},"HandleVisibility","off")
title("Pleuromamma xiphias", "FontAngle","italic"); %, "Color",soft{1})
legend(["t = 6 h", "t = 12 h"], "Location","northwest", "Box","off"); %,"TextColor", soft{1}) %, "6 h trendline", "12 h trendline"
% For the Field Data chapter.
%save("../../../2022_0214_AE2123_BC/Chemstation-S/datasets/zoopRates.mat", "pxtabs6","pxnames")


subplot(2,2,2)
% C. py plot.
cptabs = scattertabs(rI.Species == "C. pyrimidata",:);
meanmass = mean(rI.dryWeight(rI.Species == "C. pyrimidata"));
kthresh = kfun(meanmass);
%below = sum(cptabs>kthresh,1)==0;
med = median(cptabs,1,"omitnan");
[~,ia] = sort(med);
tempI = rI(rI.Species=="C. pyrimidata",:);
% idead = tempI.Notes=="DEAD";
% idead = repmat(idead',1,size(cptabs,2))';
cptabs = cptabs(:,ia);
cpnames = nicenames(~mremove); 
cpnames = cpnames(ia);
noplot = mean(cptabs,1)==0;
xCp = repelem(1:size(cpnames(~noplot),1),1,length(rI.Species(rI.Species=="C. pyrimidata")))';
bmmCp = reshape(cptabs(~noplot),length(xCp),1);
sc2 = scatter(xCp,bmmCp, 65, CP1{3},'filled',"square");
%sc2 = scatter(xCp(~idead),bmmCp(~idead), 65, CP1{3},'filled',"square");
% hold on
% sc3 = scatter(xCp(idead),bmmCp(idead), 65, CP1{3},"square");
ax = gca;
ax.YLabel.String = "pmol h^{-1} mg^{-1}";
set(ax,"XLim",[0 1+length(cpnames(~noplot))], "YScale","log");
set(ax,"XTick",1:length(cpnames(~noplot)), "XTickLabels", cpnames(~noplot),"XTickLabelRotation",45)
hold on
xL = 1.5:1:(1.5+size(cpnames(~noplot),1)-2);
h = hline(kthresh, ax.XLim);
xline(xL, "LineStyle","--", "LineWidth",0.5, "Color",[.5,.5,.5],"HandleVisibility","off")
clear tempI % idead
title("Clio pyrimidata", "FontAngle","italic")
%legend(["live","dead"], "Location","northwest")

subplot(2,2,3)
% Euph plot.
eutabs = scattertabs(rI.Species == "Euph",:);
meanmass = mean(rI.dryWeight(rI.Species == "Euph"));
kthresh = kfun(meanmass);
%below = sum(eutabs>kthresh,1)==0;
med = median(eutabs,1,"omitnan");
[~,ia] = sort(med);
tempI = rI(rI.Species=="Euph",:);
% idead = tempI.Notes=="DEAD";
% idead = repmat(idead',1,size(eutabs,2))';
eutabs = eutabs(:,ia);
noplot = mean(eutabs,1)==0;
eunames = nicenames(~mremove); 
eunames = eunames(ia);
eutabs = eutabs(:,~noplot);
eunames = eunames(~noplot);
xEu = repelem(1:size(eunames,1),1,length(rI.Species(rI.Species=="Euph")))';
bmmEu = reshape(eutabs,length(xEu),1);
sc2 = scatter(xEu,bmmEu, 65, CP1{1},'filled',"hexagram");
hold on
% sc3 = scatter(xEu(idead),bmmEu(idead), 65, CP1{1},"hexagram");
ax = gca;
ax.YLabel.String = "pmol h^{-1} mg^{-1}";
set(ax,"XLim",[0 1+length(eunames)], "YScale","log");
set(ax,"XTick",1:length(eunames), "XTickLabels", eunames,"XTickLabelRotation",45)
h = hline(kthresh, ax.XLim);
hold on
xL = 1.5:1:(1.5+size(eunames,1)-2);
xline(xL, "LineStyle","--", "LineWidth",0.5, "Color",[.5,.5,.5],"HandleVisibility","off")
clear tempI %idead
title("Hansarsia microps", "FontAngle","italic")
%legend(["\it{Hansarsia microps}","{\itStylocheiron abbreviatum} (dead)"], "Location","northwest")

subplot(2,2,4)
% Amph plot
amtabs = scattertabs(rI.Species == "Amphipod (long skinny)",:);
meanmass = mean(rI.dryWeight(rI.Species == "Amphipod (long skinny)"));
kthresh = kfun(meanmass);
%below = sum(amtabs>2.1667,1)==0;
med = median(amtabs,1,"omitnan");
[~,ia] = sort(med);
amtabs = amtabs(:,ia);
amnames = nicenames(~mremove); 
amnames = amnames(ia);
noplot = mean(amtabs,1)==0;
amnames = amnames(~noplot);
amtabs = amtabs(:,~noplot);
xAm = repelem(1:size(amnames,1),1,length(rI.Species(rI.Species=="Amphipod (long skinny)")))';
bmmAm = reshape(amtabs,length(xAm),1);
sc2 = scatter(xAm,bmmAm, 65, "k",'filled',"diamond");
hold on
ax = gca;
ax.YLabel.String = "pmol h^{-1} mg^{-1}";
set(ax,"XLim",[0 1+length(amnames)], "YScale","log");
set(ax,"XTick",1:length(amnames), "XTickLabels", amnames,"XTickLabelRotation",45)
hold on
h = hline(kthresh, ax.XLim);
xL = 1.5:1:(1.5+size(amnames,1)-2);
xline(xL, "LineStyle","--", "LineWidth",0.5, "Color",[.5,.5,.5],"HandleVisibility","off")
clear tempI idead
title("\it{Scina} spp.")

%% This section will create biomass-scaled estimates of metabolite excretion.
% It's based only on the 6 hour Pleuromamma xiphias samples and applies
% only to copepods. As it stands, it will also only do the top 8
% metabolites.

% med = median([pxtabs6],1,"omitnan");
% [~,ia] = sort(med, "descend", "MissingPlacement","last");
%FieldRates(pxtabs6(:,end-7:end), pxnames(end-7:end))

FieldRates(pxtabs6, pxnames)


%% Scatterplot; all at once.
mnames = nicenames(~mremove);
pxm = median(mtabData_pmol_mgdry_hr(~mremove,sI.Species == "PX")');
amm = median(mtabData_pmol_mgdry_hr(~mremove,sI.Species == "Amphipod (long skinny)")');
cpm = mtabData_pmol_mgdry_hr(~mremove,sI.Species == "C. pyrimidata")';
eum = median(mtabData_pmol_mgdry_hr(~mremove,sI.Species == "Euph")');
[~, ip] = sort(pxm, "descend", "MissingPlacement","last");
[~, ia] = sort(amm, "descend", "MissingPlacement","last");
[~, ic] = sort(cpm, "descend", "MissingPlacement","last");
[~, ie] = sort(eum, "descend", "MissingPlacement","last");
pnames = mnames(ip); anames = mnames(ia);
cnames = mnames(ic); enames = mnames(ie);
top10names = union(pnames(1:8),union(anames(1:8),union(cnames(1:8),enames(1:8))));
[~, ord] = ismember(top10names, mnames);

%%


scattertabs = mtabData_pmol_mgdry_hr(~mremove,t12i|t6i)';
scattertabs = scattertabs(:,ord);
rI = sI(t12i|t6i,:);
ictrl = rI.Species=="CTRL";
rI = rI(~ictrl,:);
scattertabs = scattertabs(~ictrl,:);
scattertabs(scattertabs<= 0)= NaN;
med = median(scattertabs,1,"omitnan");
[~,ia] = sort(med, "descend");
names = nicenames(~mremove);
names = names(ord);
names = names(ia);
scatter_sort = scattertabs(:,ia);
names(sum(isnan(scatter_sort),1)>9,: ) = [];
scatter_sort(:,sum(isnan(scatter_sort),1)>9) = [];
scatter_sort = flip(scatter_sort,2);
names = flip(names);
% names = names(1:15);
% scatter_sort = scatter_sort(:,1:15);



% P. xiphias plot.
pxtabs6 = scatter_sort(rI.Species == "PX"&rI.Nominal_Duration_h==6,:);
pxtabs12 = scatter_sort(rI.Species == "PX"&rI.Nominal_Duration_h==12,:);
meanmass = mean(rI.dryWeight(rI.Species == "PX"));
kthresh1 = kfun(meanmass);
xPx6 = repmat(1:size(names,1),1,sum(rI.Species=="PX"&rI.Nominal_Duration_h==6))'; %1:sum(~below)
xPx12 = repmat(1:size(names,1),1,sum(rI.Species=="PX"&rI.Nominal_Duration_h==12))'; %1:sum(~below)
bmmPx6 = [pxtabs6(1,:),pxtabs6(2,:),pxtabs6(3,:)]';
bmmPx12 = [pxtabs12(1,:),pxtabs12(2,:),pxtabs12(3,:)]';
sc1 = scatter(xPx6,bmmPx6, 65, CP1{5},"v", "LineWidth",2.5);
hold on
sc2 = scatter(xPx12,bmmPx12, 65, CP1{5},'filled',"^", "LineWidth",2.5);
xL = 1.5:1:(1.5+size(names,1)-2);
xline(xL, "LineStyle","--", "LineWidth",0.5, "Color",[.5,.5,.5],"HandleVisibility","off")


% C. py plot.
cptabs = scatter_sort(rI.Species == "C. pyrimidata",:);
meanmass = mean(rI.dryWeight(rI.Species == "C. pyrimidata"));
kthresh2 = kfun(meanmass);
tempI = rI(rI.Species=="C. pyrimidata",:);
% idead = tempI.Notes=="DEAD";
% idead = repmat(idead',1,size(cptabs,2))';
xCp = repelem(1:size(names,1),1,length(rI.Species(rI.Species=="C. pyrimidata")))';
bmmCp = reshape(cptabs,length(xCp),1);
%sc3 = scatter(xCp(~idead),bmmCp(~idead), 65, CP1{3},'filled',"square");
sc4 = scatter(xCp,bmmCp, 75, CP1{3},"filled","square", "LineWidth",2.5);
clear tempI %idead


% Euph plot.
eutabs = scatter_sort(rI.Species == "Euph",:);
meanmass = mean(rI.dryWeight(rI.Species == "Euph"));
kthresh3 = kfun(meanmass);
tempI = rI(rI.Species=="Euph",:);
% idead = tempI.Notes=="DEAD";
% idead = repmat(idead',1,size(eutabs,2))';
xEu = repelem(1:size(names,1),1,length(rI.Species(rI.Species=="Euph")))';
bmmEu = reshape(eutabs,length(xEu),1);
sc5 = scatter(xEu,bmmEu, 75, CP1{1},'filled',"hexagram", "LineWidth",2.5);
% sc6 = scatter(xEu(idead),bmmEu(idead), 65, CP1{1},"hexagram");
clear tempI idead

% Amph plot
amtabs = scatter_sort(rI.Species == "Amphipod (long skinny)",:);
meanmass = mean(rI.dryWeight(rI.Species == "Amphipod (long skinny)"));
kthresh4 = kfun(meanmass);
xAm = repelem(1:size(names,1),1,length(rI.Species(rI.Species=="Amphipod (long skinny)")))';
bmmAm = reshape(amtabs,length(xAm),1);
sc7 = scatter(xAm,bmmAm, 75, "k",'filled',"diamond", "LineWidth",2.5);
ax = gca;
set(ax,"XLim",[0 1+length(names)])
ax.YGrid = "on";
h1 = hline(kthresh1, ax.XLim);
h1.Color = CP1{5}; h1.LineWidth = 2;
h2 = hline(kthresh2, ax.XLim);
h2.Color = CP1{3}; h2.LineWidth = 2;
h2.LineStyle = ":";
h3 = hline(kthresh3, ax.XLim);
h3.Color = CP1{1}; h3.LineWidth = 2;
h4 = hline(kthresh4, ax.XLim);
h4.Color = "k"; h4.LineWidth = 2;
ax.YLabel.String = "pmol h^{-1} mg^{-1}";
set(ax, "YScale","log");
set(ax,"XTick",1:length(names), "XTickLabels", names,"XTickLabelRotation",45)
title("All Species' Top 8 Excretion Rates", "FontWeight","bold")
legend(["\it{P. xiphias}_{6h}", "\it{P. xiphias}_{12h}",...
    "\it{C. pyrimidata}",..."\it{C. pyrimidata} (dead)",...
    "\it{Hansarsia microps}",..."{\itStylocheiron abbreviatum} (dead)",...
    "\it{Scina} spp."], "Location","northwest")

if 0
    highnames = nicenames([6:6,12:12,21:24,28:28,31:31,41:41,43:43,46:46,51:53,55:56,59:60,62:65,68:69,75:75,78:78,81:81],1);
    [ihigh] = ismember(names,highnames);
    text(ax.XTick(ihigh)-0.1,0.005*ones(1,sum(ihigh)), "*", "Color","r", "FontSize", 20,"HandleVisibility","off")
end

%% An accounting of elements contained in the metabolites.

load("DOC_DON.mat")

% First: let's actually break the metabolites into elements. 
totC = mtabData_pM_Reorder_subctrl(~mremove,:).*mtabElem.C(~mremove,:);
totN = mtabData_pM_Reorder_subctrl(~mremove,:).*mtabElem.N(~mremove,:);
totO = mtabData_pM_Reorder_subctrl(~mremove,:).*mtabElem.O(~mremove,:);

elmC = totC(:,sI.Nominal_Duration_h==12);
elmN = totN(:,sI.Nominal_Duration_h==12);
elmC(isnan(elmC))= 0; elmN(isnan(elmN))= 0;
rI = sI(sI.Nominal_Duration_h==12,:);

ictrl = rI.Species=="CTRL";
rI = rI(~ictrl,:);

% meanctrlsc = mean(elmC(:,ictrl),2,"omitnan");
% meanctrlsn = mean(elmN(:,ictrl),2,"omitnan");
elmC = elmC(:,~ictrl);%-meanctrlsc;
elmN = elmN(:,~ictrl);%-meanctrlsn;
elmC(elmC<0)=0;elmN(elmN<0)=0;
itaurine = (mtabNames(~mremove)=="taurine"); iglycine = (mtabNames(~mremove)=="glycine");
ihsb = (mtabNames(~mremove)=="homoserine betaine"); iput = (mtabNames(~mremove)=="putrescine");
iarg = (mtabNames(~mremove)=="arginine"); iser = (mtabNames(~mremove)=="serine");
iala = (mtabNames(~mremove)=="alanine"); inbz = (mtabNames(~mremove)=="4-aminobenzoic acid");
Ctau = elmC(itaurine,:); Ntau = elmN(itaurine,:);
Cgly = elmC(iglycine,:); Ngly = elmN(iglycine,:);
Chsb = elmC(ihsb,:); Nhsb = elmN(ihsb,:);
Cput = elmC(iput,:); Nput = elmN(iput,:);
Carg = elmC(iarg,:); Narg = elmN(iarg,:);
% Cser = elmC(iser,:); Nser = elmN(iser,:);
Cala = elmC(iala,:); Nala = elmN(iala,:);
Cnbz = elmC(inbz,:); Nnbz = elmN(inbz,:);

elmC = sum(elmC,1);elmN = sum(elmN,1);
[G, ID] = findgroups(rI.Species);

mfunc = @(x) mean(x,1,"omitmissing");

MeanChsb = splitapply(@mean,Chsb,G'); MeanChsb = MeanChsb([4,2,1,3])./1e6;
MeanNhsb = splitapply(@mean,Nhsb,G'); MeanNhsb = MeanNhsb([4,2,1,3])./1e6;

MeanCput = splitapply(@mean,Cput,G'); MeanCput = MeanCput([4,2,1,3])./1e6;
MeanNput = splitapply(@mean,Nput,G'); MeanNput = MeanNput([4,2,1,3])./1e6;

MeanCarg = splitapply(@mean,Carg,G'); MeanCarg = MeanCarg([4,2,1,3])./1e6;
MeanNarg = splitapply(@mean,Narg,G'); MeanNarg = MeanNarg([4,2,1,3])./1e6;

% MeanCser = splitapply(@mean,Cser,G'); MeanCser = MeanCser([4,2,1,3])./1e6;
% MeanNser = splitapply(@mean,Nser,G'); MeanNser = MeanNser([4,2,1,3])./1e6;

MeanCala = splitapply(@mean,Cala,G'); MeanCala = MeanCala([4,2,1,3])./1e6;
MeanNala = splitapply(@mean,Nala,G'); MeanNala = MeanNala([4,2,1,3])./1e6;

MeanCt = splitapply(@mean,Ctau,G'); MeanCt = MeanCt([4,2,1,3])./1e6;
MeanNt = splitapply(@mean,Ntau,G'); MeanNt = MeanNt([4,2,1,3])./1e6;

MeanCg = splitapply(@mean,Cgly,G'); MeanCg = MeanCg([4,2,1,3])./1e6;
MeanNg = splitapply(@mean,Ngly,G'); MeanNg = MeanNg([4,2,1,3])./1e6;

MeanCnbz = splitapply(@mean,Cnbz,G'); MeanCnbz = MeanCnbz([4,2,1,3])./1e6;
MeanNnbz = splitapply(@mean,Nnbz,G'); MeanNnbz = MeanNnbz([4,2,1,3])./1e6;

MeanC = splitapply(@mean,elmC,G'); MeanC = MeanC([4,2,1,3])./1e6;
MeanN = splitapply(@mean,elmN,G'); MeanN = MeanN([4,2,1,3])./1e6;


mtabC = [MeanChsb;MeanCput;MeanCarg;MeanCala;MeanCt;MeanCg;MeanCnbz];% MeanCser;
mtabN = [MeanNhsb;MeanNput;MeanNarg;MeanNala;MeanNt;MeanNg;MeanNnbz];%MeanNser;
%CT = DCN{1,2:end} - DCN{1,1};
%NT = DCN{2,2:end} - DCN{2,1};
CT = DCN_ctrlsub12{1,:};
NT = DCN_ctrlsub12{2,:};
DecoyBar = CT - MeanC - sum(mtabC,1);
MeanCnomtab = MeanC - sum(mtabC,1);
DecoyBarN = NT - MeanN - sum(mtabN,1);
MeanNnomtab = MeanN - sum(mtabN,1);

load("AlbumMaps.mat","CP2")
ComboColors = flip([CP1;CP2]);

%% Updated figure to use just moles and not molar numbers
subplot(1,2,1)
b = barh(0.06.*[mtabC',MeanCnomtab',DecoyBar'], "stacked");
b(8).FaceColor = [0.2 0.2 0.2];
b(9).FaceColor = "w";
ax = gca;
order = ["Control","\it{P. xiph.}, n = 3   ", "\it{C. pyr.}, n = 1     ", "\it{Scina} spp., n = 2     ", "\it{H. microps,}  n = 2     "];
ax.YTick = [];
xlabel("DOC Relative to Control, \mumol")
set(ax, "Box", "off", "XColor", "none", "YColor","none")
ax.XDir = "reverse";
ax.XAxis.TickLabelColor = "k";
ax.XAxis.Label.Color = "k";
ax.XLim = [0 30*0.06];
ax.XGrid = "on";
percents = string(round(100.*MeanC./CT,1))+ "%";
text(b(8).YEndPoints+0.2, b(8).XEndPoints-0.05, percents, 'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
for ii = 1:length(b)-2
    b(ii).FaceColor = ComboColors{ii};
end

subplot(1,2,2)
bn = barh(0.06.*[mtabN', MeanNnomtab',DecoyBarN'], "stacked");
bn(8).FaceColor = [0.2 0.2 0.2];
bn(9).FaceColor = "w";
ax = gca;
yticklabels(order(2:5))
xlabel("TDN Relative to Control, \mumol")
set(ax, "Box", "off","XColor", "none", "YColor","none")
ax.TickLength = [0,0];
ax.TickLabelInterpreter = "tex";
ax.XAxis.TickLabelColor = "k";
ax.YAxis.TickLabelColor = "k";
ax.XAxis.Label.Color = "k";
ax.XGrid = "on";
legend(["homoserine betaine","putrescine","arginine",..."serine",...
    "alanine","taurine","glycine","4-aminobenzoic acid*",...
    "other metabolites","unaccounted"], "Location","southeast")
percentsN = string(round(100.*MeanN./NT,1))+ "%";
text(bn(8).YEndPoints+0.2, bn(8).XEndPoints-0.05, percentsN, 'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
for ii = 1:length(bn)-2
    bn(ii).FaceColor = ComboColors{ii};
end



percentCArg = 100.*MeanCarg./CT;
percentNArg = 100.*MeanNarg./NT;

