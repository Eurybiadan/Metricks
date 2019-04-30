% Copyright (C) 2019 Robert F Cooper
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
% Running this function will ask you to select the set of montage images that
% you wish to analyze. Montage images will have a single layer within a montage
% placed within an image the size of the montage canvas.
% 
% **If you only have a PSD of a montage and do not have data in this format,
% please use my PSD Layer Exporter, at 
% https://github.com/Eurybiadan/PSD_Layer_Export/releases/tag/v1.0). 
% To dump them to disk in this format.
% 
% **This function uses parfor loops, which requires the Parallel Toolbox 
%  from MATLAB.** 
%  If you don't have that toolbox, change the "parfor" on line 93 to a "for" loop.
% 
% As above, it will then prompt the user to select what the output unit
% should be. At present, the options are:
% * Microns (using millimeters^2 for density)
% * Degrees
% * Arcminutes
% 
% Once the output unit is select, it will give the user the option to pick a
% lookup table. The lookup table allows the software to analyze a folder of images
% from different subjects/timepoints/conditions. The lookup table itself **must**
% be a 3 column 'csv' file, where the **first column** is a common identifier for
% image/coordinate pairs, the **second column** is the axial length (or '24' if 
% the axial length is unknown) of the image/coordinate pairs, and the **third column** 
% is the pixels per degree of the image/coordinate pairs. Each row must contain a 
% different identifier/axial length/pixels per degree tuple.
% 
% An example common identifier could be a subject number, e.g, when working with the files
% - 1235_dateoftheyear_OD_0004.tif
% - 1235_dateoftheyear_OD_0005.tif
% 
% Common identifiers could be "1235", "1235_dateoftheyear", "1235_dateoftheyear_OD".
% If all three were placed in a LUT, then the one that matches the most (as determined
% via levenshtein distance) will be used. In this example, we would use "1235_dateoftheyear_OD".
% 
% If we had another date, say: 1235_differentdateoftheyear_OD_0005.tif, then
% _only_ the identifier "1235" would match between all images. However, say the
% two dates have different scales, then you would want to create two rows in the
% look up table for each date, with identifiers like: "1235_dateoftheyear" and 
% "1235_differentdateoftheyear".
% 
% **If you do not wish to use a lookup table, then press "cancel", and the
% software will allow you put in your own scale in UNITS/pixel.**
% 
% The software will then run, showing the spacing montage, the density montage,
% the confidence montage, and the sum map. It will also save them to disk in the
% same folder you ran from alongside a mat file that contains the results.


clear;
close all;

[fNames,thispath ]=uigetfile(fullfile(pwd,'*.tif'),'Select all files you wish to analyze from a SINGLE subject.', 'MultiSelect', 'on');

[scaling, unit, lut] = determine_scaling(pwd,fNames);

scaling = unique(scaling);
rel_scale = 1;

if length(scaling)>1
    error('The scales of the images are not the same!');
end

if ~isempty(lut) % If we didn't directly input a scale,
    % ask if we want to scale the montage to a common (smallest) scale from given LUT.
    pressedbutton = questdlg('Scale montage to the common (smallest) scale from the selected LUT?',...
                             'Scale to common size?', 'No');

    if strcmp(pressedbutton,'Yes')

        all_scales = nan(length(fNames),1);

        for i=1:length(lut{1})
            % Calculate the scale for each identifier.                                

            axiallength = lut{2}(i);
            pixelsperdegree = lut{3}(i);

            micronsperdegree = (291*axiallength)/24;

            switch unit
                case 'microns (mm density)'
                    all_scales(i) = 1 / (pixelsperdegree / micronsperdegree);
                case 'degrees'
                    all_scales(i) = 1/pixelsperdegree;
                case 'arcmin'
                    all_scales(i) = 60/pixelsperdegree;
            end
        end

        global_scale = min(all_scales); % Everything will be scaled to this.

        rel_scale = scaling./global_scale; % We will need to scale our images by this.
        scaling = global_scale; % The new scale becomes this.

    elseif strcmp(pressedbutton,'Cancel')
        return
    end
end


%% Determine the DFT distance for each image in the montage
im_spac_map = cell(length(fNames),1);
im_err_map = cell(length(fNames),1);
im_sum_map = cell(length(fNames),1);
imbox = cell(length(fNames),1);

imsize = size(imread( fullfile(thispath, fNames{1}) ));

imsize = round( imsize(1:2).*rel_scale );

fnamesplits = strsplit(fNames{1},'_');
prefix=[];
for f=1:5 % build our prefix.
    prefix = [prefix fnamesplits{f} '_'];
end

if isempty(gcp)
    myPool=parpool;
else
    myPool=gcp('nocreate');
end

parfor i=1:length(fNames)
    fNames{i}
    im = imread( fullfile(thispath, fNames{i}) );
    
    if size(im,3) >1
        im = im(:,:,1);
    end
    
    im = imresize(im, imsize);
    
    [~, im_spac_map{i}, im_err_map{i}, im_sum_map{i}, imbox{i}] = fit_fourier_spacing(im, 128);    
    
end

delete(myPool)

%%

blendedim = zeros(imsize(1:2));
blendederrim = zeros(imsize(1:2));
sum_map = zeros(imsize(1:2));

for i=1:length(imbox)
   
    thisbox = imbox{i};
    thismap = im_spac_map{i};
    thiserrmap = im_err_map{i};
    thissummap = im_sum_map{i};
    
    thismap(isnan(thismap))=0;
    thiserrmap(isnan(thiserrmap))=0;
    
    blendedim( thisbox(2):thisbox(2)+thisbox(4),...
               thisbox(1):thisbox(1)+thisbox(3) ) = blendedim( thisbox(2):thisbox(2)+thisbox(4),...
                                                               thisbox(1):thisbox(1)+thisbox(3) ) + thismap;
    
    sum_map( thisbox(2):thisbox(2)+thisbox(4),...
             thisbox(1):thisbox(1)+thisbox(3) ) = sum_map( thisbox(2):thisbox(2)+thisbox(4),...
                                                           thisbox(1):thisbox(1)+thisbox(3) ) + thissummap;

    blendederrim( thisbox(2):thisbox(2)+thisbox(4),...
               thisbox(1):thisbox(1)+thisbox(3) ) = blendederrim( thisbox(2):thisbox(2)+thisbox(4),...
                                                               thisbox(1):thisbox(1)+thisbox(3) ) + thiserrmap;                                        
end

blendedim = blendedim./blendederrim;
blendederrim = blendederrim./sum_map;


%% Display the results.
if kstest(blendederrim(~isnan(blendederrim)))
    disp('The error is not normally  distributed.')
    threshold = quantile(blendederrim(~isnan(blendederrim)),0.05)
    
else
    threshold = mean(blendederrim(:),'omitnan')-2*std(blendederrim(:),'omitnan')
    disp('The error is normally distributed.')
end

if threshold < 0.15 % Our absolute cutoff for threshold should be 0.3- that is pretty abysmal.
    threshold = 0.15
end

threshold_mask = (blendederrim>threshold);

allred = ones(256,3);
% allred(:,1) = (255:-1:0)'/255;
% allred(:,2) = (0:255)'/255;
allred(:,3) = (0:255)'/255;

% To density, assuming perfect packing
density_map = sqrt(3)./ (2*(blendedim*scaling).^2);

if strcmp(unit,'microns (mm density)')
    density_map = (1000^2).*density_map;
end
    

threshdensitymap=density_map;
rois =threshdensitymap>quantile(threshdensitymap(threshdensitymap>0),.95);

% Find our biggest cc- this will likely be the fovea.
cc = bwconncomp(rois);

numPixels = cellfun(@numel,cc.PixelIdxList);
[biggest,idx] = max(numPixels);

if length(numPixels)>=1
    numPixels = sort(numPixels,'descend')./biggest;
    % If this section is not far and away the biggest region, 
    % request that the user helps you out. Otherwise, use the biggest CC.
    if any(numPixels(2:end)>0.8)
        % To find foveal center
        fig=figure(10); clf; hold on;
        imagesc(density_map.*threshold_mask); axis image;
        title('Drag/Resize the Rectangle over the fovea.');
        h=imrect(gca,[1 1 1536 1536]);
        h.Deletable = false;
        pos= wait(h);

        roi = poly2mask([pos(1) pos(1)+pos(3) pos(1)+pos(3) pos(1) pos(1)],...
                        [pos(2) pos(2) pos(2)+pos(4) pos(2)+pos(4) pos(2)],...
                        size(density_map,1), size(density_map,2));
        close(fig)
    else
        bounding_box = regionprops(cc,'BoundingBox');
        bounding_box = bounding_box(idx).BoundingBox;
        roi = zeros(size(density_map));
        roi(cc.PixelIdxList{idx}) = 1;
    end
end

threshdensitymap=threshdensitymap.*roi;
threshdensitymap(isnan(threshdensitymap))=0;
threshdensitymap = threshdensitymap(bounding_box(2):bounding_box(2)+bounding_box(4), bounding_box(1):bounding_box(1)+bounding_box(3));
%% Find the foveal center
smoothmap = density_map;
smoothmap = smoothmap(bounding_box(2):bounding_box(2)+bounding_box(4), bounding_box(1):bounding_box(1)+bounding_box(3));
smoothmap(isnan(smoothmap)) = 0;
smoothmap = imgaussfilt(smoothmap,8);

smoothmaptheshold = quantile(smoothmap(smoothmap>0),.95);

figure(10); clf; hold on;
imagesc(smoothmap); axis image;
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
plot(foveapts(:,1),foveapts(:,2),'.');

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

fovea_coords = fovea_coords+bounding_box(1:2);


%% Make a foveal mask using active contour
smtheta_r         = linspace(0,2*pi);
smellipse_x_r     = ellipsefit.X0 + (ellipsefit.a/2)*cos( smtheta_r );
smellipse_y_r     = ellipsefit.Y0 + (ellipsefit.b/2)*sin( smtheta_r );
smrotated_ellipse = R * [smellipse_x_r;smellipse_y_r];

smfoveamask = poly2mask(smrotated_ellipse(1,:),smrotated_ellipse(2,:),size(threshdensitymap,1),size(threshdensitymap,2));
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
    
    smfoveamask = ones(size(threshdensitymap,1),size(threshdensitymap,2));
    smfoveamask(ccfovmask.PixelIdxList{mostind}) = 0;
%     figure;  imagesc(smfoveamask)
end

foveamask = ones(size(blendedim,1),size(blendedim,2));

foveamask(bounding_box(2):bounding_box(2)+bounding_box(4), bounding_box(1):bounding_box(1)+bounding_box(3)) = smfoveamask;

% foveamask = ~poly2mask(foveapts(:,1)+bounding_box(1),foveapts(:,2)+bounding_box(2),size(blendedim,1),size(blendedim,2));
% figure;  imagesc(foveamask)



%% Display and output
clear rois smoothmap
result_path = fullfile(thispath,'Results');
mkdir(result_path)


%%
save( fullfile(result_path,[prefix 'Fouriest_Result.mat']), '-v7.3' );

%% Show figures

% figure(1); imagesc(sum_map); axis image; colorbar;

blendedim = blendedim.*2/sqrt(3); % To convert to ICD spacing.
scaled_spacing = (blendedim.*scaling)-min(blendedim(:).*scaling);
scaled_spacing = 255.*(scaled_spacing./ max(scaled_spacing(:)) );

figure(2); imagesc( (blendedim.*scaling) ); axis image; colorbar;
imwrite( uint8(scaled_spacing.*threshold_mask), parula(256), fullfile(result_path, [prefix 'thresh_montage_spacing_' num2str(scaling,'%5.2f') '.tif']))
clear scaled_spacing;

scaled_error = 255*blendederrim;
[cmap, amap] = firecmap(quantile(blendederrim(blendederrim~=0), 0.01),...
                    quantile(blendederrim(blendederrim~=0), 0.25),...
                    quantile(blendederrim(blendederrim~=0), 0.05), 256);

figure(3); imagesc(blendederrim); colormap(cmap); axis image; colorbar;
imwrite( uint8(scaled_error), cmap,fullfile(result_path, [prefix 'thresh_montage_err_' num2str(scaling,'%5.2f') '.png']), 'Transparency',amap)
clear scaled_error;

imwrite( uint8(imclose(threshold_mask,ones(11)).*255), fullfile(result_path, [prefix 'thresh_montage_mask_' num2str(scaling,'%5.2f') '.tif']));

masked_density = density_map.*imclose(threshold_mask,ones(11)).*foveamask;

lower01 = quantile(masked_density(~isnan(masked_density)),0.01);
upper99 = quantile(masked_density(~isnan(masked_density)).*roi(~isnan(masked_density)),0.99);
clear masked_density

scaled_density = density_map-lower01;
scaled_density = 255.*(scaled_density./ upper99 );
scaled_density(scaled_density>255) = 255;

figure(5); imagesc( density_map.*imclose(threshold_mask,ones(11)).*foveamask ); axis image; colorbar;
caxis([lower01 upper99]);
imwrite( uint8(scaled_density.*imclose(threshold_mask,ones(11)).*foveamask),parula(256),fullfile(result_path, [prefix 'thresh_montage_density_' num2str(scaling,'%5.2f') '.tif']))
clear scaled_density;

%% Temp horizontal plot
% band_radius = 256;
% cmap = parula(256);
% 
% eccent_values = 1:length(band_avg_values);
% eccent_values = scaling*(eccent_values-fovea_coords(1));
% 
% band_avg_values = mean(density_map(fovea_coords(2)-band_radius:fovea_coords(2)+bandradius,:),'omitnan');
% 
% scaled_band_values = band_avg_values-lower01;
% scaled_band_values = 255.*(scaled_band_values./upper99);
% 
% 
% for l=1:length(band_avg_values)-1
%    
%     plot(
%     
% end




