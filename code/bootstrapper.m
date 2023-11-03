function [tr,LabelOrder,LabelsOrdered] = bootstrapper(mtabData,sI,num_boots)
% Bootstrapped Dendrogram based on metabolite data, either full dataset or
% reduced. 
% This file only uses Bray-Curtis as a distance metric, and requires you to
% have f_braycurtis.m from the Fathom toolbox in your working folder. 
% INPUT
%   mtabData: can be in any units, but dimensions are n_mtabs x n_samples
%   sI: this is an internally assumed variable, because at this point we
%       should have a sample info variable with several fields that will be
%       used to construct unique identifiers. See the Labels variable.
%   num_boots: number of times the tree will be subject to
%       reconstruction.

mn = mtabData; mn(isnan(mn))=0; % set NaNs to zero
Labels = string(sI.Species) + " " + string(round(hours(sI.duration))) +...
    "h " + string(sI.Notes) + " " + string(1:size(sI,1))'; % Sample labels.
num_mtabs = size(mn,2);
orig_mtab_dist = f_braycurtis(mn); % Generate a working distance matrix...
orig_mtab_tree = linkage(orig_mtab_dist,'average'); % and linkages
orig_mtab_tree = phytree(orig_mtab_tree, Labels);
compoundList_len = size(mn,1);

boots = cell(num_boots,1);
for n = 1:num_boots
    reorder_index = randsample(compoundList_len,compoundList_len,true);
    for i = num_mtabs:-1:1 %reverse order to preallocate memory
        bootseq(i).Header = Labels(i);
        bootseq(i).Data = mn(reorder_index,i)';
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
LabelsOrdered = LabelsOrdered(1:size(sI,1));
% This part is very specific, as it removes the whitespace and numbers that
% I placed at the end of each sample name. 
[~,LabelOrder] = ismember(LabelsOrdered, Labels);
pat = regexpPattern('\s\d+$');
LabelsOrdered = replace(LabelsOrdered, pat,'');

end