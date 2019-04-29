function [cmap, amap] = firecmap(pureredend, yellowmid, fadestart, range)

% Modal Spacing paper: 
% redcutoff: pureredend = 0.1882; % 1st percenile
% yellowmid = 0.5608; % 25th percenile
% fadestart = 0.3451; % 5th percenile

valuerange = range;
redidx = floor(valuerange*pureredend)+1; % To prevent a 0 index.
yellowidx = floor(valuerange*yellowmid)+1;

% Create our heatmap
cmap = zeros(valuerange, 3);
yellowendval = 0.8359;

gramp = 0:(yellowendval/(yellowidx-redidx)):yellowendval;
cmap(1:yellowidx,1) = 1;
cmap(redidx:yellowidx,2) = gramp;

% Create our alpha LUT
fadestartidx = floor(valuerange*fadestart)+1;

fadestartval = 1;
fadeendval = 0.3;

amap = ones(valuerange, 1);
alpharamp = fadestartval:((fadeendval-fadestartval)/...
              (yellowidx-fadestartidx)):fadeendval;
amap(fadestartidx:yellowidx) = (alpharamp);
amap(yellowidx+1:end) = 0;

figmap = repmat(amap, [1 30]);

cbar = (0:255)';
cbar = repmat(cbar,[1 30]);

% fig= figure(110); hold on;
% im=imshow(cbar,cmap);
% alpha(figmap);
% axis image; hold off; drawnow;


% imwrite(cbar, cmap, 'hot_colorbar.png', 'Transparency',amap)
end