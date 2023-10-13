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
stdectrl0 = nanstd(mtabData_pM_Reorder(:,t0ctrl_i),[],2); %/sqrt(3);
t6ctrl = mtabData_pM_Reorder(:,t6i & sInfoBig.Species=="CTRL");
meanctrl6 = nanmean(t6ctrl,2);
stdectrl6 = nanstd(t6ctrl,[],2); %/sqrt(3);
t12ctrl = mtabData_pM_Reorder(:,t12i & sInfoBig.Species=="CTRL");
meanctrl12 = nanmean(t12ctrl,2);
stdectrl12 = nanstd(t12ctrl,[],2); %/sqrt(3);
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


%% NMDS Plot of Metabolites

% AnalysisType = "Euclidean";
% AnalysisType = "Jaccard";
% AnalysisType = "BrayCurtis";
AnalysisType = "BC_BootstrapTree";
sI = sInfoBig;

switch AnalysisType
    case "Euclidean"
        mnan = mtabData_pM_Reorder; % mtabData_pmol_mgdry;
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
        title("PCA of log_2 Normalized Metabolite Data")

        % Dendrogram of the Euclidean Distance

        Labels = [string(sI.Species) + " " + string(round(hours(sI.duration))) + "h " + string(sI.Notes)];
        Z=linkage(D,'average');
        [H,T,outperm] = dendrogram(Z,size(sInfoBig,1),'Labels',Labels,'orientation','left');
        LabelOrder = flip(outperm);

    case "Jaccard"
        % Jaccard Distance with many random initial configs
        numTries = 100;
        mn = mtabData_pM_Reorder; mn(isnan(mn))=0;
        D = squareform(pdist(mn', 'jaccard'));
        opts = statset('MaxIter',10000,'Display','off','TolX',1e-5);
        Y0 = zeros(size(mn,2),2,numTries);
        stress0 = zeros(numTries,1);
        for i=1:numTries
            try
                [Y0(:,:,ii),stress0(ii)] = mdscale(D, 2, 'Start','random','Options',opts);
            catch error
                Y0(:,:,ii) = NaN;
                stress0(ii) = NaN;
            end

        end
        Y = 1e4.*mean(Y0,3, 'omitmissing');
        stress = mean(stress0, 1, 'omitmissing');

        % ml2_z = ml2;
        % ml2_z(isnan(ml2_z)|isinf(ml2_z))=0;

        sI = sInfoBig;
        sI.x = Y(:,1);
        sI.y = Y(:,2);

        px = sI.Species=="PX" & sI.Nominal_Duration_h==12;
        eu = sI.Species=="Euph";
        am = sI.Species=="Amphipod (long skinny)";
        cp = sI.Species=="C. pyrimidata";
        dead = sI.Notes=="DEAD";

        mnr = mn(sum(mn,2)>mean(sum(mn,2)),:);
        nr = mtabNames(sum(mn,2)>mean(sum(mn,2)));
        Sfc = Y(:,1:2)-mean(Y(:,1:2),1);
        mlc = [mnr',sI.dryWeight];
        mlc = mlc-mean(mlc,1,"omitnan");
        coefs = zeros(size(mlc,2),2);
        r2 = zeros(size(mlc,2),1);
        for ii=1:size(coefs,1)
            lm = fitlm(Sfc,mlc(:,ii));
            coefs(ii,:) = lm.Coefficients.Estimate(2:3)';
            r2(ii) = lm.Rsquared.Adjusted;
        end
        vecs = r2.*coefs./sqrt(sum(coefs.^2,2));
        vecnames = [nr;"Dry Weight"]; %mtabNames(r2>0.7);
        %vecs = vecs(r2>0.7,:);
        z = zeros(size(vecs,1),1);


        sc1 = scatter(sI.x(px),sI.y(px),36,stolas{1},"filled");
        hold on
        sc2 = scatter(sI.x(am),sI.y(am),36,stolas{2},"filled");
        sc3 = scatter(sI.x(eu&~dead),sI.y(eu&~dead),36,stolas{3},"filled");
        sc4 = scatter(sI.x(cp&~dead),sI.y(cp&~dead),36,stolas{4},"filled");
        sc6 = scatter(sI.x(cp&dead),sI.y(cp&dead),36,stolas{4}, "HandleVisibility","off");
        sc5 = scatter(sI.x(eu&dead),sI.y(eu&dead),36, stolas{3}, "HandleVisibility","off");
        sc7 = scatter(sI.x(t12ctrl_i), sI.y(t12ctrl_i), 36, stolas{5}, 'filled');
        q = quiver(z,z,vecs(:,1),vecs(:,2),"HandleVisibility","off");
        t = text(vecs(:,1),vecs(:,2),vecnames,"HandleVisibility","off");
        t2 = text(sI.x(cp&dead)+0.1,sI.y(cp&dead),"dead", "HandleVisibility","off");
        t3 = text(sI.x(eu&dead)+0.1,sI.y(eu&dead),"dead", "HandleVisibility","off");
        legend({'Copepod','Amphipod','Euphausiid','Pteropod','control'},"Location","best")
        xlabel("NMDS1");
        ylabel("NMDS2");
        title("NMDS of Metabolite Data, Stress = " + string(stress))

        % Dendrogram for Jaccard

        Labels = [string(sI.Species) + " " + string(round(hours(sI.duration))) + "h " + string(sI.Notes)];
        Z=linkage(D,'average');
        [H,T,outperm] = dendrogram(Z,size(sInfoBig,1),'Labels',Labels,'orientation','left');
        LabelOrder = flip(outperm);

    case "BrayCurtis"

        % Bray-Curtis with multiple random starts
        numTries = 100;
        mn = mtabData_pM_Reorder; mn(isnan(mn))=0;
        D = f_braycurtis(mn);
        opts = statset('MaxIter',10000,'Display','off','TolX',1e-5);
        Y0 = zeros(size(mn,2),2,numTries);
        stress0 = zeros(numTries,1);
        for ii=1:numTries
            [Y0(:,:,ii),stress0(ii)] = mdscale(D, 2, 'Start','random','Options',opts);
        end
        Y = 1e4.*mean(Y0,3);
        stress = mean(stress0);

        % ml2_z = ml2;
        % ml2_z(isnan(ml2_z)|isinf(ml2_z))=0;

        sI = sInfoBig;
        sI.x = Y(:,1);
        sI.y = Y(:,2);

        px = sI.Species=="PX" & sI.Nominal_Duration_h==12;
        eu = sI.Species=="Euph";
        am = sI.Species=="Amphipod (long skinny)";
        cp = sI.Species=="C. pyrimidata";
        dead = sI.Notes=="DEAD";

        mnr = mn(sum(mn,2)>mean(sum(mn,2)),:);
        nr = mtabNames(sum(mn,2)>mean(sum(mn,2)));
        Sfc = Y(:,1:2)-mean(Y(:,1:2),1);
        mlc = [mnr',sI.dryWeight];
        mlc = mlc-mean(mlc,1,"omitnan");
        coefs = zeros(size(mlc,2),2);
        r2 = zeros(size(mlc,2),1);
        for ii=1:size(coefs,1)
            lm = fitlm(Sfc,mlc(:,ii));
            coefs(ii,:) = lm.Coefficients.Estimate(2:3)';
            r2(ii) = lm.Rsquared.Adjusted;
        end
        vecs = r2.*coefs./sqrt(sum(coefs.^2,2));
        vecnames = [nr;"Dry Weight"]; %mtabNames(r2>0.7);
        %vecs = vecs(r2>0.7,:);
        z = zeros(size(vecs,1),1);

        sc1 = scatter(sI.x(px),sI.y(px),36,stolas{1},"filled");
        hold on
        sc2 = scatter(sI.x(am),sI.y(am),36,stolas{2},"filled");
        sc3 = scatter(sI.x(eu&~dead),sI.y(eu&~dead),36,stolas{3},"filled");
        sc4 = scatter(sI.x(cp&~dead),sI.y(cp&~dead),36,stolas{4},"filled");
        sc6 = scatter(sI.x(cp&dead),sI.y(cp&dead),36,stolas{4}, "HandleVisibility","off");
        sc5 = scatter(sI.x(eu&dead),sI.y(eu&dead),36, stolas{3}, "HandleVisibility","off");
        sc7 = scatter(sI.x(t12ctrl_i), sI.y(t12ctrl_i), 36, stolas{5}, 'filled');
        q = quiver(z,z,vecs(:,1),vecs(:,2),"HandleVisibility","off");
        t = text(vecs(:,1),vecs(:,2),vecnames,"HandleVisibility","off");
        t2 = text(sI.x(cp&dead)+0.1,sI.y(cp&dead),"dead", "HandleVisibility","off");
        t3 = text(sI.x(eu&dead)+0.1,sI.y(eu&dead),"dead", "HandleVisibility","off");
        legend({'Copepod','Amphipod','Euphausiid','Pteropod','control'},"Location","best")
        xlabel("NMDS1");
        ylabel("NMDS2");
        title("NMDS of Metabolite Data, Bray-Curtis Dissimilarity, Stress = " + string(stress))

        % Dendrogram for Bray-Curtis

        %distfunc = @(x,dnu) squareform(f_braycurtis(x,dnu));
        Labels = [string(sI.Species) + " " + string(round(hours(sI.duration))) + "h " + string(sI.Notes)];
        Z=linkage(D,'average');
        [H,T,outperm] = dendrogram(Z,size(sInfoBig,1),'Labels',Labels,'orientation','left');
        LabelOrder = flip(outperm);

    case "BC_BootstrapTree"
        % Bootstrapped Dendrogram
        mn = mtabData_pM_Reorder; mn(isnan(mn))=0;
        Labels = [string(sI.Species) + " " + string(round(hours(sI.duration))) + "h " + string(sI.Notes) + " " + string(1:22)'];
        num_mtabs = size(mn,2);
        orig_mtab_dist = f_braycurtis(mn);
        orig_mtab_tree = linkage(orig_mtab_dist,'average');
        orig_mtab_tree = phytree(orig_mtab_tree, Labels);
        % phytreeviewer(orig_primates_tree);

        num_boots = 1000;
        compoundList_len = size(mn,1);

        boots = cell(num_boots,1);
        for n = 1:num_boots
            reorder_index = randsample(compoundList_len,compoundList_len,true);
            for i = num_mtabs:-1:1 %reverse order to preallocate memory
                bootseq(i).Header = Labels(i);
                bootseq(i).Data = mn(:,i)';
            end
            boots{n} = bootseq;
        end

        fun = @(x,Labels) phytree(linkage(x,'average'),Labels);
        boot_trees = cell(num_boots,1);
        parpool('local');
        addAttachedFiles(gcp, "f_braycurtis.m")

        parfor (n = 1:num_boots)
            dist_tmp = f_braycurtis(vertcat(boots{n,:}.Data)');
            boot_trees{n} = fun(dist_tmp,Labels);
        end
        delete(gcp('nocreate'));

        for i = num_mtabs-1:-1:1  % for every branch, reverse order to preallocate
            branch_pointer = i + num_mtabs;
            sub_tree = subtree(orig_mtab_tree,branch_pointer);
            orig_pointers{i} = getcanonical(sub_tree);
            orig_mtabs{i} = sort(get(sub_tree,'LeafNames'));
        end

        for j = num_boots:-1:1
            for i = num_mtabs-1:-1:1  % for every branch
                branch_ptr = i + num_mtabs;
                sub_tree = subtree(boot_trees{j},branch_ptr);
                clusters_pointers{i,j} = getcanonical(sub_tree);
                clusters_mtabs{i,j} = sort(get(sub_tree,'LeafNames'));
            end
        end

        count = zeros(num_mtabs-1,1);
        for i = 1 : num_mtabs -1  % for every branch
            for j = 1 : num_boots * (num_mtabs-1)
                if isequal(orig_pointers{i},clusters_pointers{j})
                    if isequal(orig_mtabs{i},clusters_mtabs{j})
                        count(i) = count(i) + 1;
                    end
                end
            end
        end

        Pc = count ./ num_boots;   % confidence probability (Pc)

        [ptrs,dist,names] = get(orig_mtab_tree,'POINTERS','DISTANCES','NODENAMES');

        for i = 1:num_mtabs -1  % for every branch
            branch_ptr = i + num_mtabs;
            names{branch_ptr} = [names{branch_ptr} ', confidence: ' num2str(100*Pc(i)) ' %'];
        end

        tr = phytree(ptrs,dist,names);

        high_conf_branch_ptr = find(Pc > 0.9) + num_mtabs;
        view(tr, high_conf_branch_ptr)
        LabelsOrdered = get(tr,'NODENAMES');
        LabelsOrdered = LabelsOrdered(1:22);
        [~,LabelOrder] = ismember(LabelsOrdered, Labels);
        pat = regexpPattern('\s\d+$');
        LabelsOrdered = replace(LabelsOrdered, pat,'');
end

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
h = heatmap(mtabNames, LabelsOrdered+ string(1:22)', HeatMapMtabs');
h.ColorLimits = [1,10];
h.ColorScaling = "log";

colors = [stolas{1};stolas{4};stolas{5};stolas{2}];
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
    b1(1).CData = stolas{2};
    b1(2).CData = repelem(stolas{4},length(reducedNames(max(reducedRates>1,[],2))),1);
    b1(3).CData = repelem(stolas{5},length(reducedNames(max(reducedRates>1,[],2))),1);
    b1(4).CData = repelem(stolas{3},length(reducedNames(max(reducedRates>1,[],2))),1);
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

i12h = (sInfoBig.Nominal_Duration_h==12);
iDead = sInfoBig.Notes=="DEAD";

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
Px_pmol = mean(Px_pmol0,2,'omitmissing');
Px_stde_pmol = std(Px_pmol0,[],2,'omitmissing');
clear Px_pmol0

Amph_pmol0 = mtabData_pmol(:,t12i&iAmph&~iDead); Amph_pmol0(iReject_Amph) = NaN;
Amph_pmol = mean(Amph_pmol0,2,'omitmissing');
Amph_stde_pmol = std(Amph_pmol0,[],2,'omitmissing');
clear Amph_pmol0

Clio_pmol0 = mtabData_pmol(:,t12i&iClio&~iDead); Clio_pmol0(iReject_Clio) = NaN;
Clio_pmol = mean(Clio_pmol0,2,'omitmissing');
Clio_stde_pmol = std(Clio_pmol0,[],2,'omitmissing');
clear Clio_pmol0

Euph_pmol0 = mtabData_pmol(:,t12i&iEuph&~iDead); Euph_pmol0(iReject_Euph) = NaN;
Euph_pmol = mean(Euph_pmol0,2,'omitmissing');
Euph_stde_pmol = std(Euph_pmol0,[],2,'omitmissing');
clear Euph_pmol0

% Now, the time-normalized inventory.
Px_pmol_mgdry_hr0 = mtabData_pmol_mgdry_hr(:,t12i&iPx&~iDead); Px_pmol_mgdry_hr0(iReject_Px) = NaN;
Px_pmol_mgdry_hr = mean(Px_pmol_mgdry_hr0,2,'omitmissing');
Px_stde_pmol_mgdry_hr = std(Px_pmol_mgdry_hr0,[],2,'omitmissing');
clear Px_pmol_mgdry_hr0

Amph_pmol_mgdry_hr0 = mtabData_pmol_mgdry_hr(:,t12i&iAmph&~iDead); Amph_pmol_mgdry_hr0(iReject_Amph) = NaN;
Amph_pmol_mgdry_hr = mean(Amph_pmol_mgdry_hr0,2,'omitmissing');
Amph_stde_pmol_mgdry_hr = std(Amph_pmol_mgdry_hr0,[],2,'omitmissing');
clear Amph_pmol_mgdry_hr0

Clio_pmol_mgdry_hr0 = mtabData_pmol_mgdry_hr(:,t12i&iClio&~iDead); Clio_pmol_mgdry_hr0(iReject_Clio) = NaN;
Clio_pmol_mgdry_hr = mean(Clio_pmol_mgdry_hr0,2,'omitmissing');
Clio_stde_pmol_mgdry_hr = std(Clio_pmol_mgdry_hr0,[],2,'omitmissing');
clear Clio_pmol_mgdry_hr0

Euph_pmol_mgdry_hr0 = mtabData_pmol_mgdry_hr(:,t12i&iEuph&~iDead); Euph_pmol_mgdry_hr0(iReject_Euph) = NaN;
Euph_pmol_mgdry_hr = mean(Euph_pmol_mgdry_hr0,2,'omitmissing');
Euph_stde_pmol_mgdry_hr = std(Euph_pmol_mgdry_hr0,[],2,'omitmissing');
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

reducedNames = mtabNames;
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
    b1(1).CData = stolas{2};
    b1(2).CData = repelem(stolas{4},length(reducedNames(max(reducedRates>1,[],2))),1);
    b1(3).CData = repelem(stolas{5},length(reducedNames(max(reducedRates>1,[],2))),1);
    b1(4).CData = repelem(stolas{3},length(reducedNames(max(reducedRates>1,[],2))),1);
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
