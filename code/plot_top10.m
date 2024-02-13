
Px6 = t6i & iPx;
mtabData_px6 = mtabData_pmol(~mremove,Px6);
mtabData_px6(mtabData_px6<0) = NaN;
mtabData_px6m = mean(mtabData_px6,2,"omitnan");
[sortpx6m,ipx6m] = sort(mtabData_px6m, "descend", "MissingPlacement", "last");
namespx6m = nicenames(~mremove);
namespx6m = namespx6m(ipx6m);
mtabData_px6s = std(mtabData_px6,[],2,"omitnan")./sqrt(3);
sortpx6s = mtabData_px6s(ipx6m);
h = sortpx6m(1:10); n = namespx6m(1:10); er = sortpx6s(1:10);
bb = bar(1,h, "stacked");
    hold on
for ii = 1:10
    if ii==1
        yd = bb(ii).YEndPoints./2;
    else
        yd = bb(ii-1).YEndPoints + ((bb(ii).YEndPoints - bb(ii-1).YEndPoints)./2);
    end
    eb(ii) = errorbar((1.50+((-1)^ii)*0.05),yd,er(ii),er(ii), "Color",bb(ii).FaceColor);
    text(1,yd,n(ii), "HorizontalAlignment","center")
end
ylim([0,bb(10).YEndPoints+er(10)]);
xlim([0.45,1.7])
xticks([])
xlabel("6 h mean \it{P. xiphias}")
ylabel("Metabolite excreted (pmol)")
hold off