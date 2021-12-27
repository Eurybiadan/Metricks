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

minglobalbounds = round(-minglobalbounds);
%% Add all of the dft information to the montage, and combine.
avg_density = zeros(global_dimension);
% weighted_avg_spacing = zeros(global_dimension);
avg_error = zeros(global_dimension);
combined_sum_map = zeros(global_dimension,'double');

% Find the weighted averages of all of the subjects.
for f=1:length(fNames)
    disp(['Calculating Averages: ' num2str(f) ' of ' num2str(length(fNames))])
    load( fullfile(thispath, fNames{f}), 'density_map_comb', 'blendederr_comb');
    
    rowrange = round((montage_rect{f}(1,2):montage_rect{f}(3,2))+minglobalbounds(2)+1);
    colrange = round((montage_rect{f}(1,1):montage_rect{f}(3,1))+minglobalbounds(1)+1);
    
    if isempty(strfind(fNames{f}, global_eye)) % If it doesn't match our eye, flip the montage data.
        density_map_comb = fliplr(density_map_comb);
        blendederr_comb = fliplr(blendederr_comb);
    end    

    % TO DENSITY - TEMPORARY
%      blendedim = sqrt(3)./ (2*(blendedim.*scaling).^2);
%      blendedim = (1000^2).*blendedim;
    
    avg_density( rowrange, colrange) = sum(cat(3, avg_density( rowrange, colrange), density_map_comb),3,'omitnan');
    avg_error( rowrange, colrange) = sum(cat(3, avg_error( rowrange, colrange), blendederr_comb),3,'omitnan');
    combined_sum_map( rowrange, colrange) = sum(cat(3, combined_sum_map( rowrange, colrange), blendederr_comb>0),3,'omitnan');
    clear density_map_comb blendederr_comb rowrange colrange
end

avg_density = avg_density./combined_sum_map;
% avg_spacing= avg_spacing.*2/sqrt(3); % To ICD <-- only needed if row, not cell, spacing.
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
% avg_density = avg_spacing; %.*sqrt(3)/2; %<-- only needed if row, not cell, spacing.
% avg_density = (1000^2).*sqrt(3)./ (2*(avg_density).^2); %<-- for cells/mm
% avg_density = sqrt(3)./ (2*(avg_spacing).^2); %<-- for cells/deg

threshspacingmap=avg_density(minglobalbounds(2)-768:minglobalbounds(2)+768, minglobalbounds(1)-768:minglobalbounds(1)+768);

threshspacingmap(threshspacingmap<quantile(threshspacingmap(threshspacingmap>0),.85))=0;
threshspacingmap(isnan(threshspacingmap))=0;

smoothmap = imgaussfilt(threshspacingmap,8);
smoothmaptheshold= quantile(smoothmap(smoothmap>0),.85);

[clvls]=contour(smoothmap, [smoothmaptheshold smoothmaptheshold]);

[maxlvl]=max(clvls(1,:));

bloblocs = find(clvls(1,:)==maxlvl); % Sometimes we have multiple contour pieces at the same lvl.

upperclvl=[];
for b=1:length(bloblocs)
    blobval = clvls(1 ,bloblocs(b));    
    upperclvl = [upperclvl; clvls(:,bloblocs(b)+1:bloblocs(b)+clvls(2,bloblocs(b)))'];
end
convpts = convhull(upperclvl(:,1), upperclvl(:,2));
foveapts = upperclvl(convpts,:);

% foveapts = upperclvl;

% plot(foveapts(:,1),foveapts(:,2),'.');


if size(foveapts,1) > 3
    figure(10); clf; hold on;
    imagesc(smoothmap); axis image;
    plot(foveapts(:,1),foveapts(:,2),'r.'); 
    ellipsefit = fit_ellipse(foveapts(:,1),foveapts(:,2));

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


    plot( rotated_ellipse(1,:),rotated_ellipse(2,:),'r' );
    plot(fovea_coords(:,1), fovea_coords(:,2),'*');


    hold off;drawnow;

    fovea_coords = fovea_coords+minglobalbounds(1)-768;
    fovea_coords = fovea_coords+minglobalbounds(2)-768;
    
    smtheta_r         = linspace(0,2*pi);
    smellipse_x_r     = ellipsefit.X0 + (ellipsefit.a/2)*cos( smtheta_r );
    smellipse_y_r     = ellipsefit.Y0 + (ellipsefit.b/2)*sin( smtheta_r );
    smrotated_ellipse = R * [smellipse_x_r;smellipse_y_r];

    smfoveamask = poly2mask(smrotated_ellipse(1,:),smrotated_ellipse(2,:),size(threshspacingmap,1),size(threshspacingmap,2));
    smbounding = regionprops(smfoveamask,'BoundingBox');
    foveamask= activecontour(smoothmap, smfoveamask,300);

    ccfovmask = bwconncomp(foveamask);

    if ccfovmask.NumObjects >= 1
        region_info = regionprops(ccfovmask,'BoundingBox');
        mostoverlap = 0;
        mostind = 1;

        for i = 1: ccfovmask.NumObjects
            olap = rectint(region_info(i).BoundingBox, smbounding.BoundingBox);
            if olap>mostoverlap
                mostoverlap = olap;
                mostind = i;
            end
        end

        smfoveamask = ones(size(threshspacingmap,1),size(threshspacingmap,2));
        smfoveamask(ccfovmask.PixelIdxList{mostind}) = 0;
        figure;  imagesc(smfoveamask)
    end

    threshold_mask = true(size(avg_density));

    threshold_mask(minglobalbounds(2)-768:minglobalbounds(2)+768, minglobalbounds(1)-768:minglobalbounds(1)+768) = smfoveamask;

%     threshold_mask = logical(threshold_mask.*~poly2mask(splinefitx,splinefity,size(avg_spacing,1),size(avg_spacing,2)));
%     clear maxcoordsth maxcoords maxcoordsr splinefitr maxes threshspacingmap s
end
% Ensure that we only see the averages of data with more than X points going into it.
threshold_mask = logical(threshold_mask.*(combined_sum_map>5)); 

return; 
%% Plot our masked data.

spacingcaxis = quantile(avg_density(threshold_mask(:)),[0.01 0.99]);
figure(1); imagesc(avg_density.*threshold_mask); title('Combined Spacing'); axis image;
caxis(spacingcaxis); colorbar;

%% Output spacing image
rescaled_avg_spacing = (avg_density.*threshold_mask)-spacingcaxis(1);
rescaled_avg_spacing(rescaled_avg_spacing<0) = 0;
rescaled_avg_spacing = 255*(rescaled_avg_spacing./quantile(rescaled_avg_spacing(rescaled_avg_spacing~=0),0.99)); 
rescaled_avg_spacing(rescaled_avg_spacing>255) = 255;

imwrite(rescaled_avg_spacing, parula(256), [num2str(length(fNames)) 'subjects_combined_spacing.tif'])
clear rescaled_avg_spacing;

%% Output error image
errquartiles = quantile(avg_error(avg_error~=0 & ~isnan(avg_error) & ~isinf(avg_error)),[0.001 0.01 0.05]);
[firemap, amap] = firecmap(errquartiles(1), errquartiles(3), errquartiles(2),...
                           min(avg_error(avg_error~=0 & ~isnan(avg_error)& ~isinf(avg_error))), ...
                           max(avg_error(avg_error~=0 & ~isnan(avg_error)& ~isinf(avg_error))), 256);

rescaled_avg_err = avg_error.*threshold_mask;
rescaled_avg_err = 255*rescaled_avg_err; 
rescaled_avg_err(rescaled_avg_err>255) = 255;
imwrite(rescaled_avg_err, firemap, [num2str(length(fNames)) 'subjects_combined_error.png'], 'Transparency',amap)

figure(2); imagesc(rescaled_avg_err); colormap(firemap); axis image; colorbar; title('Average Error');
clear rescaled_avg_err;

%% Output density image
% avg_density = avg_density; %.*sqrt(3)/2; %<-- only needed if row, not cell, spacing.
% avg_density = (1000^2).*sqrt(3)./ (2*(avg_density).^2); %<-- for cells/mm
% avg_density = sqrt(3)./ (2*(avg_density).^2); %<-- for cells/deg

rescaled_avg_density = (avg_density.*threshold_mask)-quantile(avg_density(threshold_mask(:)),0.01);
rescaled_avg_density(rescaled_avg_density<0) = 0;
rescaled_avg_density = 255*(rescaled_avg_density./quantile(rescaled_avg_density(rescaled_avg_density~=0),0.99)); 
rescaled_avg_density(rescaled_avg_density>255) = 255;

figure(3); imagesc(avg_density.*threshold_mask); title('Combined Density'); axis image;
caxis(quantile(avg_density(threshold_mask(:)),[0.01 0.99])); colorbar;
imwrite(rescaled_avg_density, parula(256), [num2str(length(fNames)) 'subjects_combined_density.tif'])
clear rescaled_avg_density;

%% Output sum map

imwrite(uint8(255*(combined_sum_map./max(combined_sum_map(:)))), parula(256), [num2str(length(fNames)) 'subjects_sum_map.tif']);

%% Output directional plots
% avg_density = avg_density; %.*sqrt(3)/2; %<-- only needed if row, not cell, spacing.
% avg_density = (1000^2).*sqrt(3)./ (2*(avg_density).^2); %<-- for cells/mm
% avg_density = sqrt(3)./ (2*(avg_density).^2); %<-- for cells/deg

maskedspac = avg_density.*threshold_mask;
micron_position = (0:(1600/0.41))*0.41;
strip_length = length(micron_position)-1;
figure(10); clf; 
% Superior
sup_strip = mean(flipud(maskedspac(minglobalbounds(2)-strip_length:minglobalbounds(2),minglobalbounds(1)-180:minglobalbounds(1)+180)),2, 'omitnan');
plot(micron_position,sup_strip); hold on;
% Inferior
inf_strip = mean(maskedspac(minglobalbounds(2):minglobalbounds(2)+strip_length, minglobalbounds(1)-180:minglobalbounds(1)+180),2, 'omitnan');
plot(micron_position,inf_strip);
% Nasal
nasal_strip = mean(maskedspac(minglobalbounds(2)-180:minglobalbounds(2)+180,minglobalbounds(1):minglobalbounds(1)+strip_length), 1, 'omitnan');
plot(micron_position, nasal_strip);
% Temporal
temp_strip = mean(fliplr(maskedspac(minglobalbounds(2)-180:minglobalbounds(2)+180,minglobalbounds(1)-strip_length:minglobalbounds(1))), 1, 'omitnan');
plot(micron_position,temp_strip);
legend('Superior','Inferior','Nasal','Temporal')

saveas(gcf,[num2str(length(fNames)) 'subjects_directionalavg.svg']);

figure(11);
plot(micron_position, 100*(mean([sup_strip inf_strip],2)'./mean([nasal_strip; temp_strip],1)-1) );
axis([150 1600 0 50])



clear maskedspac temp_strip nasal_strip inf_strip sup_strip
return;

%% Determine average/stddev of all data.
value_map_variance = nan(global_dimension);


% For finding outliers:
figure(1); clf; hold on;
    

% Find the std deviations of all of the subjects.
for f=1:length(fNames)
    disp(['Calculating Standard Deviation: ' num2str(f) ' of ' num2str(length(fNames))])
    load( fullfile(thispath, fNames{f}), 'density_map_comb');
    
    rowrange = round((montage_rect{f}(1,2):montage_rect{f}(3,2))+minglobalbounds(2)+1);
    colrange = round((montage_rect{f}(1,1):montage_rect{f}(3,1))+minglobalbounds(1)+1);
    
    if isempty(strfind(fNames{f}, global_eye)) % If it doesn't match our eye, flip the montage data.
        density_map_comb = fliplr(density_map_comb);
    end
    
    % TO DENSITY - TEMPORARY
%     blendedim = sqrt(3)./ (2*(blendedim*scaling).^2);
%     blendedim = (1000^2).*blendedim;
    
    density_map_comb = density_map_comb.*threshold_mask( rowrange, colrange);

    thisdiff = sum( cat(3, density_map_comb, -avg_density( rowrange, colrange)) ,3).^2; % This line is good. Nans need to be carried through.
    
%     thisdiff =  sum( cat(3, (scaling.*blendedim).*(2/sqrt(3)), -avg_spacing( rowrange, colrange)) ,3).^2;

    value_map_variance( rowrange, colrange) = sum( cat(3,value_map_variance( rowrange, colrange), thisdiff), 3,'omitnan'); 
     clear thisdiff;


    clear blendedim rowrange colrange
end

% plot(micron_position,temp_strip,'k', 'LineWidth',3);
%%
 value_std_dev = value_map_variance;
 value_std_dev = value_std_dev./(combined_sum_map-1);
 value_std_dev(isinf(value_std_dev)) =0;
 value_std_dev = sqrt(value_std_dev);
 maskedspac = value_std_dev.*threshold_mask;

figure(3); imagesc(value_std_dev.*threshold_mask); title('Combined Spacing Std dev');
saveas(gcf,[num2str(length(fNames)) 'subjects_combined_stddev.svg']);

figure(20); clf; hold on;
% Superior
plot(micron_position/291,mean(flipud(maskedspac(minglobalbounds(2)-strip_length:minglobalbounds(2),minglobalbounds(1)-180:minglobalbounds(1)+180)),2,  'omitnan')); hold on;
% Inferior
plot(micron_position/291,mean(maskedspac(minglobalbounds(2):minglobalbounds(2)+strip_length, minglobalbounds(1)-180:minglobalbounds(1)+180),2, 'omitnan'));
% Nasal
nasal_strip = mean(maskedspac(minglobalbounds(2)-180:minglobalbounds(2)+180,minglobalbounds(1):minglobalbounds(1)+strip_length), 1,  'omitnan');
plot(micron_position/291, nasal_strip);
% Temporal
temp_strip = mean(fliplr(maskedspac(minglobalbounds(2)-180:minglobalbounds(2)+180,minglobalbounds(1)-strip_length:minglobalbounds(1))), 1,  'omitnan');
plot(micron_position/291,temp_strip);
legend('Superior','Inferior','Nasal','Temporal')
saveas(gcf,[num2str(length(fNames)) 'subjects_directionalstddev.svg']);
clear maskedspac temp_strip nasal_strip
return;


%% Determine weighted average/stddev of all data.
weighted_avg_spacing = weighted_avg_spacing./weighted_avg_error;
weighted_avg_spacing_std_dev = zeros(global_dimension);
weighted_avg_error_std_dev = zeros(global_dimension);


% Find the weighted std deviations of all of the subjects.
for f=1:length(fNames)
    load( fullfile(thispath, fNames{f}), 'blendedim', 'blendederrim');
    
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
