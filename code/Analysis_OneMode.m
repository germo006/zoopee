% Noah Germolus 25 August 2022
% This is a quick analysis script for the Zoop data. 

clear; clc;  
close all

load("../datasets/zoopee_pM_OneMode.mat")

setDefaultFigs

load("AlbumMaps.mat", "stolas")

%% Now on to the actual data.
% Notes from lab notebook:
%   

Info = '../datasets/mtab_Noah_Zoop2_BC_071423.xlsx';
Info = readtable(Info);
AmyInfo = '../datasets/2023_May_small_boat_metabolomics_filled.xlsx';
AmyInfo = readtable(AmyInfo, 'ReadVariableNames', true);

bigInfo = join(Info(~isnan(Info.bottleNumAmy),:), AmyInfo(~isnan(AmyInfo.BTL_ID),:), 'LeftKeys', {'bottleNumAmy'}, 'RightKeys', {'BTL_ID'});

[sInfoBig, iBig] = innerjoin(sInfo, bigInfo, ...
    'LeftKeys', 'FileName_neg', 'RightKeys', 'File_Name');

% Reorder metabolite data 
mtabData_pM_Reorder = mtabData_pM(:,iBig);

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


sInfoBig.duration = sInfoBig.TimeStop - sInfoBig.TimeStart;
t0ctrl_i = find(sInfoBig.Nominal_Duration_h == 0);
t6ctrl_i = find(sInfoBig.Nominal_Duration_h == 6 & sInfoBig.Species == "CTRL");
t12ctrl_i = find(sInfoBig.Nominal_Duration_h == 12 & sInfoBig.Species == "CTRL");
meanctrl0 = mean(mtabData_pM_Reorder(:,t0ctrl_i),2);
meanctrl6 = mean(mtabData_pM_Reorder(:,t6ctrl_i),2);
meanctrl12 = mean(mtabData_pM_Reorder(:,t12ctrl_i),2);
stdvctrl = std(mtabData_pM_Reorder(:,t0ctrl_i),[],2);
mtabData_pM_Reorder_subctrl = mtabData_pM_Reorder;
mtabData_pM_Reorder_subctrl(:,sInfoBig.Nominal_Duration_h==0) = mtabData_pM_Reorder_subctrl(:,sInfoBig.Nominal_Duration_h==0) - meanctrl0;
mtabData_pM_Reorder_subctrl(:,sInfoBig.Nominal_Duration_h==6) = mtabData_pM_Reorder_subctrl(:,sInfoBig.Nominal_Duration_h==6) - meanctrl6;
mtabData_pM_Reorder_subctrl(:,sInfoBig.Nominal_Duration_h==12) = mtabData_pM_Reorder_subctrl(:,sInfoBig.Nominal_Duration_h==12) - meanctrl12;
mtabData_pmol = (mtabData_pM_Reorder_subctrl' .* (sInfoBig.Volume_mL_./1000))';
mtabData_pmol_mgdry = (mtabData_pmol' ./ (sInfoBig.dryWeight))';
mtabData_pmol_mgdry_hr = (mtabData_pmol_mgdry' ./ hours(sInfoBig.duration))';

%% Formatting the data for the figs

outdir = '../figs/zoopfigs1';
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

% I need to make these data into something easy to graph. 
% Begin with the concentration data straight up. 
% Controls (3 at t0 and tf, two singletons in the middle)

% Useful indices for time points and animals.
t0i = sInfoBig.Nominal_Duration_h == 0;
t6i = sInfoBig.Nominal_Duration_h == 6;
t12i = sInfoBig.Nominal_Duration_h==12;
iClio = sInfoBig.Species == "C. pyrimidata";
iPx = sInfoBig.Species == "PX";
iAmph = sInfoBig.Species == "Amphipod (long skinny)";
iEuph = sInfoBig.Species == "Euph";
iAnimal = iClio | iPx | iAmph | iEuph;


t0ctrl_i = find(t0i);
meanctrl0 = nanmean(mtabData_pM_Reorder(:,t0ctrl_i),2);
stdectrl0 = nanstd(mtabData_pM_Reorder(:,t0ctrl_i),[],2)/sqrt(3);
t6ctrl = mtabData_pM_Reorder(:,t6i & sInfoBig.Species=="CTRL");
meanctrl6 = nanmean(t6ctrl,2);
stdectrl6 = nanstd(t6ctrl,[],2)/sqrt(3);
t12ctrl = mtabData_pM_Reorder(:,t12i & sInfoBig.Species=="CTRL");
meanctrl12 = nanmean(t12ctrl,2);
stdectrl12 = nanstd(t12ctrl,[],2)/sqrt(3);
ctrls = [meanctrl0, meanctrl6, meanctrl12];
ste_ctrls = [stdectrl0, stdectrl6, stdectrl12];

% Times! Yeah, I should probably put some bars on that. 
meantimes = hours([mean(sInfoBig.duration(t0i)),...
    mean(sInfoBig.duration(t6i)),...
    mean(sInfoBig.duration(t12i))]);
rangetimes = hours([range(sInfoBig.duration(t0i)),...
    range(sInfoBig.duration(t6i)),...
    range(sInfoBig.duration(t12i))]);

% Now, the straight concentrations. 
% For this type of data, I'll deal with the P.xiph data first since it has
% the most complicated form. I won't need to repeat the time vector again. 
Pxtimes = meantimes(2:end);
Pxtimer = rangetimes(2:end);

Px_pM = [nanmean(mtabData_pM_Reorder(:,t6i&iPx),2),...
    nanmean(mtabData_pM_Reorder(:,t12i&iPx),2)];
Px_stde_pM = [nanstd(mtabData_pM_Reorder(:,t6i&iPx),[],2),...
    nanstd(mtabData_pM_Reorder(:,t12i&iPx),[],2)]./sqrt(3);

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
    nanstd(mtabData_pM_Reorder_subctrl(:,t12i&iPx),[],2)]./sqrt(3);

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
    nanstd(mtabData_pmol(:,t12i&iPx),[],2)]./sqrt(3);

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
    nanstd(mtabData_pmol_mgdry_hr(:,t12i&iPx),[],2)]./sqrt(3);

Amph_pmol_mgdry_hr = nanmean(mtabData_pmol_mgdry_hr(:,iAmph),2);
Amph_stde_pmol_mgdry_hr = nanstd(mtabData_pmol_mgdry_hr(:,iAmph),[],2);

Clio_pmol_mgdry_hr = nanmean(mtabData_pmol_mgdry_hr(:,iClio),2);
Clio_stde_pmol_mgdry_hr = nanstd(mtabData_pmol_mgdry_hr(:,iClio),[],2);

Euph_pmol_mgdry_hr = nanmean(mtabData_pmol_mgdry_hr(:,iEuph),2);
Euph_stde_pmol_mgdry_hr = nanstd(mtabData_pmol_mgdry_hr(:,iEuph),[],2);


%% NMDS Plot of Metabolites

mnan = mtabData_pmol_mgdry;
mnan(mnan<=0) = NaN;
minm = min(min(mnan,[],2,"omitmissing"),[],1,"omitmissing");
c = round(log2(minm)); d = 2^c;
mnan(isnan(mnan))=0;
ml2 = log2(mnan + d) - c;
D = squareform(pdist(ml2',@naneucdist));
[Y,eig] = cmdscale(D);
pv = eig./sum(eig);

% ml2_z = ml2;
% ml2_z(isnan(ml2_z)|isinf(ml2_z))=0;

Sfc = Y(:,1:2)-mean(Y(:,1:2),1);
mlc = ml2'-mean(ml2',1,"omitnan");
coefs = zeros(size(mlc,2),2);
r2 = zeros(size(mlc,2),1);
for ii=1:length(coefs)
    lm = fitlm(Sfc,mlc(:,ii));
    coefs(ii,:) = lm.Coefficients.Estimate(2:3)';
    r2(ii) = lm.Rsquared.Adjusted;
end
vecs = coefs./sqrt(sum(coefs.^2,2));
vecnames = mtabNames(r2>0.7);
vecs = vecs(r2>0.7,:);
z = zeros(size(vecs,1),1);

sI = sInfoBig;
sI.x = Y(:,1);
sI.y = Y(:,2);

px = sI.Species=="PX" & sI.Nominal_Duration_h==12;
eu = sI.Species=="Euph";
am = sI.Species=="Amphipod (long skinny)";
cp = sI.Species=="C. pyrimidata";
dead = sI.Notes=="DEAD";


sc1 = scatter(sI.x(px),sI.y(px),36,stolas{1},"filled");
hold on
sc2 = scatter(sI.x(am),sI.y(am),36,stolas{2},"filled");
sc3 = scatter(sI.x(eu&~dead),sI.y(eu&~dead),36,stolas{3},"filled");
sc4 = scatter(sI.x(cp&~dead),sI.y(cp&~dead),36,stolas{4},"filled");
sc6 = scatter(sI.x(cp&dead),sI.y(cp&dead),36,stolas{4}, "HandleVisibility","off");
sc5 = scatter(sI.x(eu&dead),sI.y(eu&dead),36, stolas{3}, "HandleVisibility","off");
q = quiver(z,z,100*vecs(:,1),100*vecs(:,2),"HandleVisibility","off");
t = text(100*vecs(:,1),100*vecs(:,2),vecnames,"HandleVisibility","off");
t2 = text(sI.x(cp&dead)+1,sI.y(cp&dead),"dead", "HandleVisibility","off");
t3 = text(sI.x(eu&dead)+1,sI.y(eu&dead),"dead", "HandleVisibility","off");
legend({'Copepod','Amphipod','Euphausiid','Pteropod'},"Location","best")
xlabel(["PC1 ("+string(round(100*pv(1),2))+"%)"]); 
ylabel(["PC2 ("+string(round(100*pv(2),2))+"%)"]);
title("NMDS of log_2 Normalized Metabolite Data")

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
            'Color',stolas{1},'LineWidth',1.5, 'LineStyle', 'none');
        eb2 = errorbar(ax1, Pxtimes, Px_pM(i,:),...
            Px_stde_pM(i,:),Px_stde_pM(i,:),...
            'Color',stolas{3},'LineWidth',1.5, 'LineStyle', 'none');
        eb3 = errorbar(ax1, meantimes(3), Amph_pM(i,:),...
            Amph_stde_pM(i,:),Amph_stde_pM(i,:),...
            'Color',stolas{4},'LineWidth',1.5, 'LineStyle', 'none');
        eb4 = errorbar(ax1, meantimes(3), Clio_pM(i,:),...
            Clio_stde_pM(i,:),Clio_stde_pM(i,:),...
            'Color',stolas{5},'LineWidth',1.5, 'LineStyle', 'none');
        eb5 = errorbar(ax1, meantimes(3), Euph_pM(i,:),...
            Euph_stde_pM(i,:),Euph_stde_pM(i,:),...
            'Color',stolas{2},'LineWidth',1.5, 'LineStyle', 'none');
        ylim([0,max([1,Px_pM(i,:),Amph_pM(i,:),Clio_pM(i,:),Euph_pM(i,:)])])
      
        
        % The legend.
        legend(ax1, {'LowLim','HighLim','ctrl','{\it P. xiphias}','{\it Amphipoda spp.}',...
            '{\it C. pyrimidata}', '{\it Euphasiid spp.}'}, 'Orientation', 'horizontal',...
            'Position', [0.1 0.48 0.82 0.04])
        
        
        % The control subtraction.
        eb6 = errorbar(ax2, meantimes, ctrls(i,:)-meanctrl0(i),...
            ste_ctrls(i,:),ste_ctrls(i,:),...
            'Color',stolas{1},'LineWidth',1.5, 'LineStyle', 'none');
        hold(ax2)
        eb7 = errorbar(ax2, Pxtimes, Px_pM_subctrl(i,:),...
            Px_stde_pM(i,:),Px_stde_pM(i,:),...
            'Color',stolas{3},'LineWidth',1.5, 'LineStyle', 'none');
        eb8 = errorbar(ax2, meantimes(3), Amph_pM_subctrl(i,:),...
            Amph_stde_pM(i,:),Amph_stde_pM(i,:),...
            'Color',stolas{4},'LineWidth',1.5, 'LineStyle', 'none');
        eb9 = errorbar(ax2, meantimes(3), Clio_pM_subctrl(i,:),...
            Clio_stde_pM(i,:),Clio_stde_pM(i,:),...
            'Color',stolas{5},'LineWidth',1.5, 'LineStyle', 'none');
        eb10 = errorbar(ax2, meantimes(3), Euph_pM_subctrl(i,:),...
            Euph_stde_pM(i,:),Euph_stde_pM(i,:),...
            'Color',stolas{2},'LineWidth',1.5, 'LineStyle', 'none'); 
        ylim([0,max([1,Px_pM_subctrl(i,:),Amph_pM_subctrl(i,:),Clio_pM_subctrl(i,:),Euph_pM_subctrl(i,:)])])

        % The inventory.
        eb11 = errorbar(ax3, meantimes, (ctrls(i,:)-meanctrl0(i)).*0.06,...
            ste_ctrls(i,:).*0.06,ste_ctrls(i,:).*0.06,...
            'Color',stolas{1},'LineWidth',1.5, 'LineStyle', 'none');
        hold(ax3)
        eb12 = errorbar(ax3, Pxtimes, Px_pmol(i,:),...
            Px_stde_pmol(i,:),Px_stde_pmol(i,:),...
            'Color',stolas{3},'LineWidth',1.5, 'LineStyle', 'none');
        eb13 = errorbar(ax3, meantimes(3), Amph_pmol(i,:),...
            Amph_stde_pmol(i,:),Amph_stde_pmol(i,:),...
            'Color',stolas{4},'LineWidth',1.5, 'LineStyle', 'none');
        eb14 = errorbar(ax3, meantimes(3), Clio_pmol(i,:),...
            Clio_stde_pmol(i,:),Clio_stde_pmol(i,:),...
            'Color',stolas{5},'LineWidth',1.5, 'LineStyle', 'none');
        eb15 = errorbar(ax3, meantimes(3), Euph_pmol(i,:),...
            Euph_stde_pmol(i,:),Euph_stde_pmol(i,:),...
            'Color',stolas{2},'LineWidth',1.5, 'LineStyle', 'none');
        ylim([0,max([1,Px_pmol(i,:),Amph_pmol(i,:),Clio_pmol(i,:),Euph_pmol(i,:)])])
        
        % The mass-normalized inventory.
        eb16 = errorbar(ax4, meantimes, ctrls(i,:),...
            ste_ctrls(i,:),ste_ctrls(i,:),...
            'Color',stolas{1},'LineWidth',1.5, 'LineStyle', 'none', 'Visible', 'off');
        hold(ax4)
        eb17 = errorbar(ax4, Pxtimes, Px_pmol_mgdry_hr(i,:),...
            Px_stde_pmol_mgdry_hr(i,:),Px_stde_pmol_mgdry_hr(i,:),...
            'Color',stolas{3},'LineWidth',1.5, 'LineStyle', 'none');
        eb18 = errorbar(ax4, meantimes(3), Amph_pmol_mgdry_hr(i,:),...
            Amph_stde_pmol_mgdry_hr(i,:),Amph_stde_pmol_mgdry_hr(i,:),...
            'Color',stolas{4},'LineWidth',1.5, 'LineStyle', 'none');
        eb19 = errorbar(ax4, meantimes(3), Clio_pmol_mgdry_hr(i,:),...
            Clio_stde_pmol_mgdry_hr(i,:),Clio_stde_pmol_mgdry_hr(i,:),...
            'Color',stolas{5},'LineWidth',1.5, 'LineStyle', 'none');
        eb20 = errorbar(ax4, meantimes(3), Euph_pmol_mgdry_hr(i,:),...
            Euph_stde_pmol_mgdry_hr(i,:),Euph_stde_pmol_mgdry_hr(i,:),...
            'Color',stolas{2},'LineWidth',1.5, 'LineStyle', 'none');
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

reducedNames = mtabNames;%(goodStuff);
elem = mtabElem; 
%elem = elem(goodStuff,:);

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
    f1 = figure('Name',"AllEstimates",'Color','none','Position',...
        [50 30 800 1000],'Units','inches');
    
    b1 = bar(categorical(reducedNames(max(reducedRates>1,[],2))),...
        reducedRates(reducedRates>1), 'FaceColor','flat');
    % b1(1).XData = categorical(reducedNames(max(reducedRates>1,[],2)));
    % b1(2).XData = categorical(reducedNames(max(reducedRates>1,[],2)));
    % b1(3).XData = categorical(reducedNames(max(reducedRates>1,[],2)));
    % b1(4).XData = categorical(reducedNames(max(reducedRates>1,[],2)));
    b1(1).CData = stolas{2};
    b1(2).CData = repelem(stolas{4},length(reducedNames(max(reducedRates>1,[],2))),1);
    b1(3).CData = repelem(stolas{5},length(reducedNames(max(reducedRates>1,[],2))),1);
    b1(4).CData = repelem(stolas{3},length(reducedNames(max(reducedRates>1,[],2))),1);
    set(gca, 'XTickLabel', reducedNames(max(reducedRates>1,[],2)), 'XTick', 1:length(reducedNames(max(reducedRates>1,[],2))),...
        'XTickLabelRotation', 90, 'YScale', 'log', 'TickLength', [0 0], ...
        'YGrid', 'on', 'YLim', [1, max(max(reducedRates(reducedRates>1)))+10]);
    ylabel('pmol mg^{-1} h^{-1}', 'Interpreter', 'tex')
    legend(gca, {'{\it P. xiphias}','{\it Amphipoda spp.}',...
        '{\it C. cuspidata}', '{\it Euphasiid spp.}'}, 'Orientation', 'vertical','Location',...
        'northwest', 'FontSize', 8)
    saveas(f1, outdir + "/"+"EstimatedRates.pdf")
end
%% Analysis of Rate Data. 

rateEst = table();
rateEst.mtabName = mtabNames;
i12h = (sInfoBig.Nominal_Duration_h==12);
rateEst.Px = nanmean(mtabData_pmol_mgdry_hr(:,i12h & iPx),2);
rateEst.Clio = nanmean(mtabData_pmol_mgdry_hr(:,i12h & iClio),2);
rateEst.Amph = nanmean(mtabData_pmol_mgdry_hr(:,i12h & iAmph),2);
rateEst.Euph = nanmean(mtabData_pmol_mgdry_hr(:,i12h & iEuph),2);
rateEst.PxFlag = nanmean(mtabData_pM_Reorder(:,i12h & iPx),2)>MaxStd_pM;
rateEst.ClioFlag = nanmean(mtabData_pM_Reorder(:,i12h & iClio),2)>MaxStd_pM;
rateEst.AmphFlag = nanmean(mtabData_pM_Reorder(:,i12h & iAmph),2)>MaxStd_pM;
rateEst.EuphFlag = nanmean(mtabData_pM_Reorder(:,i12h & iEuph),2)>MaxStd_pM;
rateEst.ErrPx = Px_stde_pmol_mgdry_hr(:,2);
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


writetable(rateEst, '../datasets/rateEstimates_24Aug.xlsx')
