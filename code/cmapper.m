function [cmap] = cmapper(colormat,n)
%CMAPPER creates a gradient (n x 3) of three-element color triplets. 
%   the colormat input is simply any number of color triplets that you want
%   to interpolate, stacked in an m x 3 configuration. 
%   The exact size of the output isn't exact--the function rounds based on
%   the number of interpolated colors and the number of points expected.
numBreaks = size(colormat,1);
ptsPerSection = round(n/(numBreaks-1));

cmap = zeros(n,3);
ind = 0;
for ii = 1:numBreaks-1
cmap(ind*ptsPerSection+1:(1+ind)*ptsPerSection,:) =...
    [linspace(colormat(ii,1),colormat(ii+1,1),ptsPerSection)',...
    linspace(colormat(ii,2),colormat(ii+1,2),ptsPerSection)'...
    linspace(colormat(ii,3),colormat(ii+1,3),ptsPerSection)'];
ind = ind + 1;
end

cmap(sum(cmap,2)==0,:) = [];