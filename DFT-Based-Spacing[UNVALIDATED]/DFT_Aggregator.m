% Robert F Cooper 2018-11-07
%
%
clear;
close all force;
[fNames,thispath ]=uigetfile(fullfile(pwd,'*.mat'),'Select all montages you wish to combine.', 'MultiSelect', 'on');

foveal_coords = nan(length(fNames), 2);
montage_size = nan(length(fNames), 2);
montage_rect = cell(length(fNames), 1);

for f=1:length(fNames)
    load( fullfile(thispath, fNames{f}), 'fovea_coords', 'imsize');
    
    foveal_coords(f,:) = fovea_coords;
    montage_size(f,:) = imsize(1:2);
    
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
weighted_avg_spacing = zeros(global_dimension);
weighted_avg_error = zeros(global_dimension);
combined_sum_map = zeros(global_dimension);

% Find the weighted averages of all of the subjects.
for f=1:length(fNames)
    load( fullfile(thispath, fNames{f}), 'fovea_coords', 'imbox', 'blendedim',...
                                         'blendederrim', 'density_map','sum_map');
    
    rowrange = round((montage_rect{f}(1,2):montage_rect{f}(3,2))-minglobalbounds(2)+1);
    colrange = round((montage_rect{f}(1,1):montage_rect{f}(3,1))-minglobalbounds(1)+1);
    
    weighted_avg_spacing( rowrange, colrange) = sum(cat(3, weighted_avg_spacing( rowrange, colrange), blendedim),3,'omitnan');
    weighted_avg_error( rowrange, colrange) = sum(cat(3, weighted_avg_error( rowrange, colrange), blendederrim),3,'omitnan');
    combined_sum_map( rowrange, colrange) = sum(cat(3, combined_sum_map( rowrange, colrange), blendederrim>0),3,'omitnan');
end

weighted_avg_spacing = weighted_avg_spacing./weighted_avg_error;


weighted_avg_spacing_std_dev = zeros(global_dimension);
weighted_avg_error_std_dev = zeros(global_dimension);


% Find the weighted std deviations of all of the subjects.
for f=1:length(fNames)
    load( fullfile(thispath, fNames{f}), 'fovea_coords', 'imbox', 'blendedim',...
                                         'blendederrim', 'density_map','sum_map');
    
    rowrange = round((montage_rect{f}(1,2):montage_rect{f}(3,2))-minglobalbounds(2)+1);
    colrange = round((montage_rect{f}(1,1):montage_rect{f}(3,1))-minglobalbounds(1)+1);
    
    % Weighted std deviation equation found here
    % https://stats.stackexchange.com/questions/6534/how-do-i-calculate-a-weighted-standard-deviation-in-excel
    weighted_avg_spacing_std_dev( rowrange, colrange) = blendederrim.*(sum(cat(3, blendedim, -weighted_avg_spacing( rowrange, colrange)),3,'omitnan').^2);

end

weighted_std_dev_err_sum = ((combined_sum_map-1)./combined_sum_map).*weighted_avg_error;
weighted_avg_spacing_std_dev  = sqrt(weighted_avg_spacing_std_dev./weighted_std_dev_err_sum);
weighted_avg_spacing_std_dev(isinf(weighted_avg_spacing_std_dev)) =0 ;
weighted_avg_error = weighted_avg_error./combined_sum_map;

% Display the results.
if kstest(weighted_avg_error(~isnan(weighted_avg_error)))
    disp('The error is not normally  distributed.')
    threshold = quantile(weighted_avg_error(~isnan(weighted_avg_error)),0.025)
    
else
    threshold = mean(weighted_avg_error(:),'omitnan')-2*std(weighted_avg_error(:),'omitnan')
    disp('The error is normally distributed.')
end

threshold_mask = (weighted_avg_error>threshold);

figure(1); imagesc(weighted_avg_spacing.*threshold_mask); title('Combined Spacing');

figure(2); imagesc(weighted_avg_spacing_std_dev.*threshold_mask); title('Combined Spacing Weighted Std dev');
figure(3); imagesc(weighted_avg_error.*threshold_mask);
