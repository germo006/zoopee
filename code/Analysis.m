% Noah Germolus 28 Nov 2023
% This is meant to be the slimmed-down analysis script which puts out just
% the results and figures that will be necessary. 

clear; clc;  
close all

% Loads data and metadata.
loadstuff

outdir = '../figs/zoopfigs1';
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

% Information for constructing a color gradient later on. 
colors = [CP1{2};CP1{1};CP1{3}];
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

clear rmNames mtab_pM_notctrl mtab_pM_ctrl

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
LiveDead = sI.Notes(~iallctrl); LiveDead = string(LiveDead(LabelOrder));
TimePointInfo = sI.Nominal_Duration_h(~iallctrl); 
TimePointInfo = string(TimePointInfo(LabelOrder));
hxdata = nicenames(~mremove);
hydata = LabelsOrdered+ " " + TimePointInfo +" h " + LiveDead+ " " +...
    string(1:size(LiveDead,1))';
hxdata(ihrm,:)=[];
HeatMapMtabs(ihrm,:)=[];
% Make the heatmap.
h = heatmap(hydata, hxdata, HeatMapMtabs);
% h.ColorLimits = [1,10];
h.ColorScaling = "scaledcolumns";
h.Colormap = flip(cmapper(colors,npts));
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


%% Scatterplot; by species and only the t12 samples

scattertabs = mtabData_pmol_mgdry_hr(~mremove,t12i|t6i)';
scattertabs(isnan(scattertabs))= 0;
rI = sI(t12i|t6i,:);
ictrl = rI.Species=="CTRL";
rI = rI(~ictrl,:);
scattertabs = scattertabs(~ictrl,:);
scattertabs(scattertabs<0)=0;
% ia12 = rI.Nominal_Duration_h==12&(rI.Species=="PX"|rI.Species=="C. pyrimidata"|rI.Species=="Euph"|rI.Species=="Amphipod (long skinny)");
scattertabs = scattertabs./rI.dryWeight./hours(rI.duration);
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
[~,ia] = sort(med);
% pxtabs6 = pxtabs6(:,ia(~below));
% pxtabs12 = pxtabs12(:,ia(~below));
pxnames = nicenames(~mremove); 
pxnames = pxnames(ia); %pxnames(ia(~below));
xPx6 = repmat(1:size(pxnames,1),1,sum(rI.Species=="PX"&rI.Nominal_Duration_h==6))'; %1:sum(~below)
xPx12 = repmat(1:size(pxnames,1),1,sum(rI.Species=="PX"&rI.Nominal_Duration_h==12))'; %1:sum(~below)
bmmPx6 = reshape(pxtabs6(:,ia),length(xPx6),1);
bmmPx12 = reshape(pxtabs12(:,ia),length(xPx12),1);
sc1 = scatter(xPx6,bmmPx6, 40, CP1{2},"v");
hold on
sc2 = scatter(xPx12,bmmPx12, 40, CP1{2},'filled',"^");
ax = gca;
ax.YLabel.String = "pmol h^{-1} mg^{-1}";
set(ax,"XLim",[0 1+length(pxnames)], "YScale","log");
set(ax,"XTick",1:length(pxnames), "XTickLabels", pxnames,"XTickLabelRotation",45)
h = hline(kthresh, ax.XLim);
xL = 1.5:1:(1.5+size(pxnames,1)-2);
xline(xL, "LineStyle","--", "LineWidth",0.5, "Color",[.5,.5,.5],"HandleVisibility","off")
title("Pleuromamma xiphias", "FontAngle","italic")
legend(["t = 6 h", "t = 12 h"], "Location","northwest")

subplot(2,2,2)
% C. py plot.
cptabs = scattertabs(rI.Species == "C. pyrimidata",:);
meanmass = mean(rI.dryWeight(rI.Species == "C. pyrimidata"));
kthresh = kfun(meanmass);
below = sum(cptabs>kthresh,1)==0;
med = median(cptabs,1,"omitnan");
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
h = hline(kthresh, ax.XLim);
xline(xL, "LineStyle","--", "LineWidth",0.5, "Color",[.5,.5,.5],"HandleVisibility","off")
clear tempI idead
title("Clio pyrimidata", "FontAngle","italic")
legend(["live","dead"], "Location","northwest")

subplot(2,2,3)
% Euph plot.
eutabs = scattertabs(rI.Species == "Euph",:);
meanmass = mean(rI.dryWeight(rI.Species == "Euph"));
kthresh = kfun(meanmass);
below = sum(eutabs>kthresh,1)==0;
med = median(eutabs,1,"omitnan");
[~,ia] = sort(med);
tempI = rI(rI.Species=="Euph",:);
idead = tempI.Notes=="DEAD";
idead = repmat(idead',1,sum(~below))';
eutabs = eutabs(:,ia(~below));
eunames = nicenames(~mremove); 
eunames = eunames(ia(~below));
xEu = repelem(1:sum(~below),1,length(rI.Species(rI.Species=="Euph")))';
bmmEu = reshape(eutabs,length(xEu),1);
sc2 = scatter(xEu(~idead),bmmEu(~idead), 40, CP1{1},'filled',"hexagram");
hold on
sc3 = scatter(xEu(idead),bmmEu(idead), 40, CP1{1},"hexagram");
ax = gca;
ax.YLabel.String = "pmol h^{-1} mg^{-1}";

set(ax,"XLim",[0 1+length(eunames)], "YScale","log");
set(ax,"XTick",1:length(eunames), "XTickLabels", eunames,"XTickLabelRotation",45)
h = hline(kthresh, ax.XLim);
hold on
xL = 1.5:1:(1.5+sum(~below)-2);
xline(xL, "LineStyle","--", "LineWidth",0.5, "Color",[.5,.5,.5],"HandleVisibility","off")
clear tempI idead
title("Euphausiids")
legend(["\it{Hansarsia microps}","{\itStylocheiron abbreviatum} (dead)"], "Location","northwest")

subplot(2,2,4)
% Amph plot
amtabs = scattertabs(rI.Species == "Amphipod (long skinny)",:);
meanmass = mean(rI.dryWeight(rI.Species == "Amphipod (long skinny)"));
kthresh = kfun(meanmass);
below = sum(amtabs>2.1667,1)==0;
med = median(amtabs,1,"omitnan");
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
h = hline(kthresh, ax.XLim);
xL = 1.5:1:(1.5+sum(~below)-2);
xline(xL, "LineStyle","--", "LineWidth",0.5, "Color",[.5,.5,.5],"HandleVisibility","off")
clear tempI idead
title("Scina spp.", "FontAngle","italic")

%% An accounting of elements contained in the metabolites.

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






