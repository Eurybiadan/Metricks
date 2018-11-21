% Robert F Cooper 2018-11-07
%
%
clear;
close all force;
[fNames,thispath ]=uigetfile(fullfile(pwd,'*.mat'),'Select all montages you wish to combine.', 'MultiSelect', 'on');

if ~all(iscell(fNames))
    fNames = {fNames};
end

for f=1:length(fNames)
    load( fullfile(thispath, fNames{f}), 'scaling');
    scales(f) = scaling;   
end
scaling = unique(scales);

if length(scaling)>1
    error('The montage scales must be the same!');
end

global_eye = 'OD'; % Everything will be flipped to this.

foveal_coords = nan(length(fNames), 2);
montage_size = nan(length(fNames), 2);
montage_rect = cell(length(fNames), 1);

for f=1:length(fNames)
    load( fullfile(thispath, fNames{f}), 'fovea_coords', 'imsize');
    
    montage_size(f,:) = imsize(1:2);
    foveal_coords(f,:) = fovea_coords;
    
    if isempty(strfind(fNames{f}, global_eye)) % If it doesn't match our eye, flip the fovea coordinates in the x direction.
        foveal_coords(f,1) = montage_size(f,2)-foveal_coords(f,1);
    end
    
    montage_rect{f} = [               1                  1; % TLC
                       montage_size(f,2)                 1; % TRC
                       montage_size(f,2) montage_size(f,1); % BRC
                                      1  montage_size(f,1); % BLC
                                      1                  1];% TLC 
    
    % Shift all montage rectangles to a null point
    montage_rect{f} = montage_rect{f}-foveal_coords(f,:);    
end

% Find the bounding rectangle
minglobalbounds = min(cell2mat(montage_rect));
maxglobalbounds = max(cell2mat(montage_rect));

global_dimension = fliplr(ceil(maxglobalbounds-minglobalbounds))+1; % Flipped so height/width (row/col), not width/height (x/y)


%% Add all of the dft information to the montage, and combine.
avg_spacing = zeros(global_dimension);
% weighted_avg_spacing = zeros(global_dimension);
avg_error = zeros(global_dimension);
combined_sum_map = zeros(global_dimension,'double');

% Find the weighted averages of all of the subjects.
for f=1:length(fNames)
    load( fullfile(thispath, fNames{f}), 'blendedim', 'blendederrim', 'sum_map');
    
    rowrange = round((montage_rect{f}(1,2):montage_rect{f}(3,2))-minglobalbounds(2)+1);
    colrange = round((montage_rect{f}(1,1):montage_rect{f}(3,1))-minglobalbounds(1)+1);
    
    if isempty(strfind(fNames{f}, global_eye)) % If it doesn't match our eye, flip the montage data.
        blendedim = fliplr(blendedim);
        blendederrim = fliplr(blendederrim);
    end
    
    avg_spacing( rowrange, colrange) = sum(cat(3, avg_spacing( rowrange, colrange), blendedim),3,'omitnan');
%     weighted_avg_spacing( rowrange, colrange) = sum(cat(3, weighted_avg_spacing( rowrange, colrange), blendedim.*blendederrim),3,'omitnan');
    avg_error( rowrange, colrange) = sum(cat(3, avg_error( rowrange, colrange), blendederrim),3,'omitnan');
    combined_sum_map( rowrange, colrange) = sum(cat(3, combined_sum_map( rowrange, colrange), blendederrim>0),3,'omitnan');
    clear blendedim blendederrim rowrange colrange
end

avg_spacing = scaling.*avg_spacing./combined_sum_map;
avg_spacing= avg_spacing.*2/sqrt(3); % To ICD
avg_error = avg_error./combined_sum_map;

% Display the results.
if kstest(avg_error(~isnan(avg_error)))
    disp('The error is not normally distributed.')
    threshold = quantile(avg_error(~isnan(avg_error)),0.05)
    
else
    threshold = mean(avg_error(:),'omitnan')-2*std(avg_error(:),'omitnan')
    disp('The error is normally distributed.')
end

threshold_mask = (avg_error>threshold) & (combined_sum_map>10);

figure(1); imagesc(avg_spacing.*threshold_mask); title('Combined Spacing');
figure(2); imagesc(avg_error.*threshold_mask); colormap(flipud(jet(256))); axis image; colorbar; title('Average Error');
return;

%% Determine average/stddev of all data.
spacing_std_dev = zeros(global_dimension);


% Find the std deviations of all of the subjects.
for f=1:length(fNames)
    load( fullfile(thispath, fNames{f}), 'blendedim','sum_map');
    
    rowrange = round((montage_rect{f}(1,2):montage_rect{f}(3,2))-minglobalbounds(2)+1);
    colrange = round((montage_rect{f}(1,1):montage_rect{f}(3,1))-minglobalbounds(1)+1);
    
    if isempty(strfind(fNames{f}, global_eye)) % If it doesn't match our eye, flip the montage data.
        blendedim = fliplr(blendedim);
    end
    
    spacing_std_dev( rowrange, colrange) = (sum(cat(3, scaling.*blendedim, -avg_spacing( rowrange, colrange)),3,'omitnan').^2);
    clear blendedim sum_map rowrange colrange
end

spacing_std_dev = sqrt(spacing_std_dev./(combined_sum_map-1));



% figure(3); imagesc(spacing_std_dev.*threshold_mask); title('Combined Spacing Std dev');

return;

%% Determine weighted average/stddev of all data.
weighted_avg_spacing = weighted_avg_spacing./weighted_avg_error;
weighted_avg_spacing_std_dev = zeros(global_dimension);
weighted_avg_error_std_dev = zeros(global_dimension);


% Find the weighted std deviations of all of the subjects.
for f=1:length(fNames)
    load( fullfile(thispath, fNames{f}), 'blendedim', 'blendederrim', 'sum_map');
    
    rowrange = round((montage_rect{f}(1,2):montage_rect{f}(3,2))-minglobalbounds(2)+1);
    colrange = round((montage_rect{f}(1,1):montage_rect{f}(3,1))-minglobalbounds(1)+1);
    
    if isempty(strfind(fNames{f}, global_eye)) % If it doesn't match our eye, flip the montage data.
        blendedim = fliplr(blendedim);
        blendederrim = fliplr(blendederrim);
    end
    
    % Weighted std deviation equation found here
    % https://stats.stackexchange.com/questions/6534/how-do-i-calculate-a-weighted-standard-deviation-in-excel
    weighted_avg_spacing_std_dev( rowrange, colrange) = blendederrim.*(sum(cat(3, blendedim, -weighted_avg_spacing( rowrange, colrange)),3,'omitnan').^2);

end

weighted_std_dev_err_sum = ((combined_sum_map-1)./combined_sum_map).*weighted_avg_error;
weighted_avg_spacing_std_dev  = sqrt(weighted_avg_spacing_std_dev./weighted_std_dev_err_sum);
weighted_avg_spacing_std_dev(isinf(weighted_avg_spacing_std_dev)) =0 ;



figure(10); imagesc(weighted_avg_spacing.*threshold_mask); title('Combined Weighted Spacing');
figure(11); imagesc(weighted_avg_spacing_std_dev.*threshold_mask); title('Combined Weighted Spacing Weighted Std dev');
