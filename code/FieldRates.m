%% Noah Germolus 12 Feb 2024
% This is a code that does what Amy's allometric respiration code for
% MOCNESS data does, but worse. 

night = readtable("../datasets/M9_test_Classified_FOR_MATLAB.csv");
day = readtable("../datasets/M8_test_Classified_FOR_MATLAB.csv");

night.Taxa = string(night.Taxa); day.Taxa = string(day.Taxa);

night(contains(night.Taxa, "not-living"), :) = [];
day(contains(day.Taxa, "not-living"), :) = [];

temp = zeros(2,1);
percents = table(temp, temp, temp, temp, ...
    'VariableNames', {'Copepoda','Amphipoda','Euphausiacea','Thecosomata'},...
    'RowNames',{'Day','Night'});
clear temp
perfunc = @(Taxon, Data) sum(contains(Data.Taxa,Taxon))/height(Data);

for ii = 1:length(percents.Properties.VariableNames)
    percents{1,ii} = perfunc(percents.Properties.VariableNames(ii),day);
    percents{2,ii} = perfunc(percents.Properties.VariableNames(ii),night);
end

% Okay, copepods are 70 % of all these MOCNESS finds, so let's just use
% copepods. 

night(~contains(night.Taxa, "Copepoda"), :) = [];
day(~contains(day.Taxa, "Copepoda"), :) = [];

night(contains(night.Taxa, "dead"), :) = [];
day(contains(day.Taxa, "dead"), :) = [];

Copepee_avg_day = day.DW*median(pxtabs6, 1, "omitmissing");
Copepee_std_day = day.DW*std(pxtabs6, 1, "omitmissing");
Copepee_avg_night = night.DW*median(pxtabs6, 1, "omitmissing");
Copepee_std_night = night.DW*std(pxtabs6, 1, "omitmissing");

sizefrac_day = findgroups(day.net,day.fraction);
sizefrac_night = findgroups(night.net,night.fraction);

Copepee_avg_day_fracs = splitapply(@sum,Copepee_avg_day,sizefrac_day);
Copepee_avg_night_fracs = splitapply(@sum,Copepee_avg_night,sizefrac_night);

Copepee_std_day_fracs = splitapply(@sum,Copepee_std_day,sizefrac_day);
Copepee_std_night_fracs = splitapply(@sum,Copepee_std_night,sizefrac_night);

clear Copepee_avg_day Copepee_std_day Copepee_avg_night Copepee_std_night

reducedDay = day(:,["moc", "net", "fraction", "D_N", "split", "Min_depth","Max_depth","hdif", "Tow_Vol"]);
reducedNight = night(:,["moc", "net", "fraction", "D_N", "split", "Min_depth","Max_depth","hdif", "Tow_Vol"]);

reducedDay = unique(reducedDay,"stable");
reducedNight = unique(reducedNight, "stable");

Copepee_avg_day_fracs = Copepee_avg_day_fracs./reducedDay.split;
Copepee_avg_night_fracs = Copepee_avg_night_fracs./reducedNight.split;
Copepee_std_day_fracs = Copepee_std_day_fracs./reducedDay.split;
Copepee_std_night_fracs = Copepee_std_night_fracs./reducedNight.split;

reducedDay.xsecV = reducedDay.hdif.*reducedDay.Tow_Vol;
reducedNight.xsecV = reducedNight.hdif.*reducedNight.Tow_Vol;

dayNets = reducedDay(reducedDay.fraction=="d3",:);
nightNets = reducedNight(reducedNight.fraction=="d3",:);

nets = [dayNets;nightNets]; nets.split = [];

addgrps_night = findgroups(reducedNight.net);
addgrps_day = findgroups(reducedDay.net);

cpeeday = splitapply(@sum,Copepee_avg_day_fracs,addgrps_day);
cpeeday_std = splitapply(@sum,Copepee_std_day_fracs,addgrps_day);
cpeenight = splitapply(@sum,Copepee_avg_night_fracs,addgrps_night);
cpeenight_std = splitapply(@sum,Copepee_std_night_fracs,addgrps_night);

cpeeavg = 12.*[cpeeday;cpeenight]./nets.Tow_Vol;
cpeestd = 12.*[cpeeday_std;cpeenight_std]./nets.Tow_Vol;

id = sum(cpeeavg,1) ==0;
colNames = pxnames(~id);
cpeeavg(:,id) = [];
cpeestd(:,id) = [];

%%
% [a ia] = sort(median(cpeeavg), "descend");
%%
%mtab = 7; 
mtab = find(colNames == "putrescine");
f = figure("Position",[100 100 2000 1000]);
ax1 = axes(f, "Position",[0.1300,0.1100,0.3695,0.8150]);
ax2 = axes(f,"Position", [0.5005, 0.1100, 0.3695, 0.8150]);
b1 = barh(ax1,...
    mean([nets.Min_depth(nets.D_N=="d"),nets.Max_depth(nets.D_N=="d")],2),...
    cpeeavg(nets.D_N=="d",mtab));
b2 = barh(ax2,...
    mean([nets.Min_depth(nets.D_N=="n"),nets.Max_depth(nets.D_N=="n")],2),...
    cpeeavg(nets.D_N=="n",mtab));
ax1.YDir = "reverse"; ax2.YDir = "reverse";
ax1.XDir = "reverse"; ax2.YAxisLocation = "right";
ax1.Box = "off"; ax2.Box = "off";
ax1.TickLength = [0 0];
ax2.TickLength = [0 0];
ax1.Color = CP1{2};
ax2.Color = CP1{3};
b1.FaceColor = CP1{4};
b2.FaceColor = CP1{1};
b1.EdgeColor = [1.00,0.93,0.81];
b2.EdgeColor = [0.58,0.72,0.81];
b1.LineWidth = 2;b2.LineWidth = 2;
% ax1.YAxis.Color = "none";
% ax2.YAxis.Color = "none";
ax1.YLim = ax2.YLim;
xl = max([ax1.XLim;ax2.XLim]);
ax1.XLim = xl; ax2.XLim = xl;
ax1.YLabel.String = "Avg Net Depth, m";
ax2.XLabel.String = "12-hour concentration increase, pM";
ax1.XTick(1) = [];
title(ax1, colNames(mtab))
ax1.FontSize = 14; ax2.FontSize = 14;
defont = "Corbel"; 
ax1.FontName = defont; ax2.FontName = defont;
ax1.Title.FontSize = 24;
t1 = text(ax1, 7*ax1.XLim(2)/8, 0,  "day");
t1.FontName = defont; t1.Color = b1.FaceColor; t1.FontSize = 24;
t1.HorizontalAlignment = "left";t1.FontWeight = "bold";
t2 = text(ax2, 7*ax2.XLim(2)/8, 0,  "night");
t2.FontName = defont; t2.Color = b1.EdgeColor; t2.FontSize = 24;
t2.HorizontalAlignment = "right"; t2.FontWeight = "bold";