function cdata = splitcolors(C, cat)
cat = categorical(cat);
G = findgroups(cat);
cdata = zeros(size(G,1),3);
for i = 1:size(G,1)
    cdata(i,:) = C{G(i)};
end
end