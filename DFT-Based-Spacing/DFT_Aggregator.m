% Copyright (C) 2019 Robert F Cooper, created 2018-11-07
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
    disp(['Calculating Foveal Location: ' num2str(f) ' of ' num2str(length(fNames))])
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

minglobalbounds = -minglobalbounds;
%% Add all of the dft information to the montage, and combine.
avg_spacing = zeros(global_dimension);
% weighted_avg_spacing = zeros(global_dimension);
avg_error = zeros(global_dimension);
combined_sum_map = zeros(global_dimension,'double');

% Find the weighted averages of all of the subjects.
for f=1:length(fNames)
    disp(['Calculating Averages: ' num2str(f) ' of ' num2str(length(fNames))])
    load( fullfile(thispath, fNames{f}), 'blendedim', 'blendederrim', 'sum_map');
    
    rowrange = round((montage_rect{f}(1,2):montage_rect{f}(3,2))+minglobalbounds(2)+1);
    colrange = round((montage_rect{f}(1,1):montage_rect{f}(3,1))+minglobalbounds(1)+1);
    
    if isempty(strfind(fNames{f}, global_eye)) % If it doesn't match our eye, flip the montage data.
        blendedim = fliplr(blendedim);
        blendederrim = fliplr(blendederrim);
    end    

    % TO DENSITY - TEMPORARY
%     blendedim = sqrt(3)./ (2*(blendedim.*scaling).^2);
%     blendedim = (1000^2).*blendedim;
    
    avg_spacing( rowrange, colrange) = sum(cat(3, avg_spacing( rowrange, colrange), blendedim),3,'omitnan');
    avg_error( rowrange, colrange) = sum(cat(3, avg_error( rowrange, colrange), blendederrim),3,'omitnan');
    combined_sum_map( rowrange, colrange) = sum(cat(3, combined_sum_map( rowrange, colrange), blendederrim>0),3,'omitnan');
    clear blendedim blendederrim rowrange colrange
end

% avg_spacing = avg_spacing./combined_sum_map; %TO DENSITY-TEMPORARY
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

threshold_mask = (avg_error>threshold); % & (combined_sum_map>10);

%% To find foveal mask
threshspacingmap=avg_spacing(minglobalbounds(2)-768:minglobalbounds(2)+768, minglobalbounds(1)-768:minglobalbounds(1)+768);
quantile(threshspacingmap(threshspacingmap>0),.15)
threshspacingmap(threshspacingmap>quantile(threshspacingmap(threshspacingmap>0),.15))=0;
threshspacingmap(isnan(threshspacingmap))=0;

maxes = imregionalmin(threshspacingmap);
% maxes = imregionalmax(threshspacingmap);%TO DENSITY-TEMPORARY

s = regionprops(maxes,'centroid');
maxcoords = zeros(length(s),2);
for i=1:length(s)
    maxcoords(i,:)  = s(i).Centroid;
end

if size(maxcoords,1) > 3
    figure(10); clf; hold on;
    imagesc(threshspacingmap); axis image;
    plot(maxcoords(:,1),maxcoords(:,2),'b.'); 
    ellipsefit = fit_ellipse(maxcoords(:,1),maxcoords(:,2));
    [maxcoordsth, maxcoordsr] = cart2pol(maxcoords(:,1)-ellipsefit.X0_in,maxcoords(:,2)-ellipsefit.Y0_in);
    f=fit(maxcoordsth,maxcoordsr,'smoothingspline','SmoothingParam',0.9992513623689557);


    maxcoordsth = sort(maxcoordsth);
    splinefitr = f(maxcoordsth);
    [splinefitx, splinefity]= pol2cart(maxcoordsth, splinefitr);
    splinefitx = splinefitx + ellipsefit.X0_in;
    splinefity = splinefity + ellipsefit.Y0_in;

    fovea_coords = [ellipsefit.X0_in ellipsefit.Y0_in];


    % rotation matrix to rotate the axes with respect to an angle phi
    cos_phi = cos( ellipsefit.phi );
    sin_phi = sin( ellipsefit.phi );
    R = [ cos_phi sin_phi; -sin_phi cos_phi ];

    % the ellipse
    theta_r         = linspace(0,2*pi);
    ellipse_x_r     = ellipsefit.X0 + ellipsefit.a*cos( theta_r );
    ellipse_y_r     = ellipsefit.Y0 + ellipsefit.b*sin( theta_r );
    rotated_ellipse = R * [ellipse_x_r;ellipse_y_r];


    plot( splinefitx, splinefity,'g.')
    plot( rotated_ellipse(1,:),rotated_ellipse(2,:),'r' );

    hold off;drawnow;

    splinefitx = splinefitx+minglobalbounds(1)-768;
    splinefity = splinefity+minglobalbounds(2)-768;
    threshold_mask = logical(threshold_mask.*~poly2mask(splinefitx,splinefity,size(avg_spacing,1),size(avg_spacing,2)));
    clear maxcoordsth maxcoords maxcoordsr splinefitr maxes threshspacingmap s
end

%% Plot our masked data.

spacingcaxis = quantile(avg_spacing(threshold_mask(:)),[0.01 0.99]);
figure(1); imagesc(avg_spacing.*threshold_mask); title('Combined Spacing'); axis image;
caxis(spacingcaxis); colorbar;

%% Output spacing image
rescaled_avg_spacing = (avg_spacing.*threshold_mask)-spacingcaxis(1);
rescaled_avg_spacing(rescaled_avg_spacing<0) = 0;
rescaled_avg_spacing = 255*(rescaled_avg_spacing./quantile(rescaled_avg_spacing(rescaled_avg_spacing~=0),0.99)); 
rescaled_avg_spacing(rescaled_avg_spacing>255) = 255;

imwrite(rescaled_avg_spacing, parula(256), [num2str(length(fNames)) 'subjects_combined_spacing.tif'])
clear rescaled_avg_spacing;

%% Output error image
errquartiles = quantile(avg_error(avg_error~=0),[0.01 0.05 0.25]);
[firemap, amap] = firecmap(errquartiles(1), errquartiles(3), errquartiles(2), 256);

rescaled_avg_err = avg_error;
rescaled_avg_err = 255*avg_error; 
rescaled_avg_err(rescaled_avg_err>255) = 255;
imwrite(rescaled_avg_err, firemap, [num2str(length(fNames)) 'subjects_combined_error.png'], 'Transparency',amap)

figure(2); imagesc(rescaled_avg_err); colormap(firemap); axis image; colorbar; title('Average Error');
clear rescaled_avg_err;

%% Output density image
avg_density = avg_spacing.*sqrt(3)/2;
avg_density = (1000^2).*sqrt(3)./ (2*(avg_density).^2);

rescaled_avg_density = (avg_density.*threshold_mask)-quantile(avg_density(threshold_mask(:)),0.01);
rescaled_avg_density(rescaled_avg_density<0) = 0;
rescaled_avg_density = 255*(rescaled_avg_density./quantile(rescaled_avg_density(rescaled_avg_density~=0),0.99)); 
rescaled_avg_density(rescaled_avg_density>255) = 255;

figure(3); imagesc(avg_density.*threshold_mask); title('Combined Density'); axis image;
caxis(quantile(avg_density(threshold_mask(:)),[0.01 0.99])); colorbar;
imwrite(rescaled_avg_density, parula(256), [num2str(length(fNames)) 'subjects_combined_density.tif'])
clear avg_density rescaled_avg_density;

%% Output sum map

imwrite(uint8(255*(combined_sum_map./max(combined_sum_map(:)))), parula(256), [num2str(length(fNames)) 'subjects_sum_map.tif']);

%% Output directional plots
maskedspac = avg_spacing.*threshold_mask;
micron_position = (0:(1600/0.41))*0.41;
strip_length = length(micron_position)-1;
figure(10); clf; 
% Superior
sup_strip = mean(flipud(maskedspac(minglobalbounds(2)-strip_length:minglobalbounds(2),minglobalbounds(1)-180:minglobalbounds(1)+180)),2);
plot(micron_position,sup_strip); hold on;
% Inferior
inf_strip = mean(maskedspac(minglobalbounds(2):minglobalbounds(2)+strip_length, minglobalbounds(1)-180:minglobalbounds(1)+180),2);
plot(micron_position,inf_strip);
% Nasal
nasal_strip = mean(maskedspac(minglobalbounds(2)-180:minglobalbounds(2)+180,minglobalbounds(1):minglobalbounds(1)+strip_length), 1);
plot(micron_position, nasal_strip);
% Temporal
temp_strip = mean(fliplr(maskedspac(minglobalbounds(2)-180:minglobalbounds(2)+180,minglobalbounds(1)-strip_length:minglobalbounds(1))), 1);
plot(micron_position,temp_strip);
legend('Superior','Inferior','Nasal','Temporal')

figure(11);
plot(micron_position, 100*(mean([sup_strip inf_strip],2)'./mean([nasal_strip; temp_strip],1)-1) );
axis([150 1600 0 50])

clear maskedspac temp_strip nasal_strip inf_strip sup_strip
return;

%% Determine average/stddev of all data.
spacing_std_dev = zeros(global_dimension);


% Find the std deviations of all of the subjects.
for f=1:length(fNames)
    disp(['Calculating Standard Deviation: ' num2str(f) ' of ' num2str(length(fNames))])
    load( fullfile(thispath, fNames{f}), 'blendedim','sum_map');
    
    sum_map=sum_map>0;
    
    rowrange = round((montage_rect{f}(1,2):montage_rect{f}(3,2))+minglobalbounds(2)+1);
    colrange = round((montage_rect{f}(1,1):montage_rect{f}(3,1))+minglobalbounds(1)+1);
    
    if isempty(strfind(fNames{f}, global_eye)) % If it doesn't match our eye, flip the montage data.
        blendedim = fliplr(blendedim);
        sum_map = fliplr(sum_map);
    end
    
    % TO DENSITY - TEMPORARY
%     blendedim = sqrt(3)./ (2*(blendedim*scaling).^2);
%     blendedim = (1000^2).*blendedim;
%     spacing_std_dev( rowrange, colrange) = spacing_std_dev( rowrange, colrange)+ ...
%                                        sum( cat(3, (blendedim.*sum_map), ...
%                                                     -avg_spacing( rowrange, colrange).*sum_map) ,3,'omitnan').^2;                                                

    spacing_std_dev( rowrange, colrange) = spacing_std_dev( rowrange, colrange)+ ...
                                           sum( cat(3, (scaling.*blendedim.*sum_map).*(2/sqrt(3)), ...
                                                        -avg_spacing( rowrange, colrange).*sum_map) ,3,'omitnan').^2;

    clear blendedim sum_map rowrange colrange
end

spacing_std_dev = sqrt(spacing_std_dev./(combined_sum_map-1));


figure(3); imagesc(spacing_std_dev.*threshold_mask); title('Combined Spacing Std dev');

maskedspac = spacing_std_dev.*threshold_mask;
figure(20); clf; 
% Superior
plot(micron_position,mean(flipud(maskedspac(minglobalbounds(2)-strip_length:minglobalbounds(2),minglobalbounds(1)-180:minglobalbounds(1)+180)),2)); hold on;
% Inferior
plot(micron_position,mean(maskedspac(minglobalbounds(2):minglobalbounds(2)+strip_length, minglobalbounds(1)-180:minglobalbounds(1)+180),2));
% Nasal
nasal_strip = mean(maskedspac(minglobalbounds(2)-180:minglobalbounds(2)+180,minglobalbounds(1):minglobalbounds(1)+strip_length), 1);
plot(micron_position, nasal_strip);
% Temporal
temp_strip = mean(fliplr(maskedspac(minglobalbounds(2)-180:minglobalbounds(2)+180,minglobalbounds(1)-strip_length:minglobalbounds(1))), 1);
plot(micron_position,temp_strip);
clear maskedspac temp_strip nasal_strip
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


lower01 = quantile(weighted_avg_spacing(~isnan(weighted_avg_spacing)),0.01);
upper99 = quantile(weighted_avg_spacing(~isnan(weighted_avg_spacing)),0.99);

figure(10); imagesc(weighted_avg_spacing.*threshold_mask); title('Combined Weighted Spacing');
caxis([lower01 upper99]);
figure(11); imagesc(weighted_avg_spacing_std_dev.*threshold_mask); title('Combined Weighted Spacing Weighted Std dev');
