function [ output_args ] = large_vessel_mask( image )
% large_vessel_mask(image)
%
% Masks out large vessels from AO images using a entropy-based approach.
%
% Created 2017-12-12


nhood_size = 31;
nhood_radius = floor(nhood_size/2);

nhood = strel('disk', nhood_radius);

% Entropy and gaussian filter the image. Sections without much entropy will
% be darker (more uniform) areas.

logim = log10(double(image+1));
logim = logim - min(logim(:));
logim = 255*logim ./ max(logim(:));

entropyim = entropyfilt(uint8(logim), nhood.getnhood);
entgausim = imgaussfilt(entropyim, nhood_radius);
entgausim = entgausim - min(entgausim(:));
entgausim = entgausim ./ max(entgausim(:));

threshim = entgausim;

% Determine our threshold.
figure(1);
gaushist = histogram(threshim,'Normalization','probability');

percent_threshold = .0075;
overthresh = find(gaushist.Values>=percent_threshold);
ourthresh = gaushist.BinEdges(overthresh(1))

% ourthresh = (2*sqrt( mean( (threshim(:)-.87).^2)))

vesselmask = threshim>ourthresh;

figure(2); imagesc(logim); colormap gray; axis image;
figure(3); imagesc(entgausim); colormap gray; axis image;
figure(4); imagesc(vesselmask); colormap gray; axis image;


end