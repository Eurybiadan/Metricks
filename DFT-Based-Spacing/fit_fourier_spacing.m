function [avg_pixel_spac, interped_spac_map, interped_conf_map, sum_map, imbox ] = fit_fourier_spacing(test_image, roi_size, supersampling, row_or_cell, roi_step)
% FUNCTION [avg_pixel_spac, interped_spac_map, interped_conf_map, sum_map, imbox ] = fit_fourier_spacing(test_image, roi_size, supersampling, row_or_cell, roi_step)
% 
% #### Inputs:
% 
% - **test_image**: The image that will be analyzed. The only requirement is that it is a 2d, grayscale (1 channel) image.
% - **roi_size**: The side length (in pixels) of a sliding roi window- The roi will march along the image you've provided at a rate of 1/roi_step the size of the ROI, creating a "map" of spacing of the image.
% - **supersampling**: If "true", then each roi will be super-sampled in accordance with: [Bernstein et al.](https://arxiv.org/pdf/1401.2636.pdf) before calculating the DFT-derived spacing.
% - **row_or_cell**: The range of angles from the polar DFT that will be used to calculate the DFT-derived spacing. If "row", then it will be the upper and lower 90 degrees of the DFT. If "cell", it will be the left and right 90 degrees.
% - **roi_step**: The step size of our sliding window. If undefined, it will default to 1/4 of the roi_size.
% 
% #### Outputs:
% 
% - **avg_pixel_spac**: The average spacing of the image.
% - **interped_spac_map**: The spacing map of the input image (in pixel spacing).
% - **interped_conf_map**: The confidence map of the input image.
% - **sum_map**: The map corresponding to the amount of ROI overlap across the output map.
% - **imbox**: The bounding region of valid (nonzero, NaN, or Inf) pixels.
%
%
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

if ~exist('test_image','var') || isempty(test_image)
    [filename, pathname] = uigetfile('*.tif', 'Pick an image to segment');

    test_image = imread( fullfile(pathname, filename) );
    
end
% tic;
if size(test_image,3) >1
    test_image = test_image(:,:,1);
end
im_size = size(test_image);

if ~exist('roi_size','var') 
    roi_size = im_size;
end

if ~exist('supersampling','var')
    supersampling = false;
end

if ~exist('row_or_cell','var')
    row_or_cell = 'cell';
end

if ~exist('roi_step','var')
    roi_step = floor(roi_size/4);
end
interped_spac_map=[];

imcomps = bwconncomp( imclose(test_image>0,ones(5)) );
imbox = regionprops(imcomps, 'BoundingBox');


boxsizes = zeros(size(imbox,1),1);
for i=1:size(imbox,1)
    boxsizes(i)= imbox(i).BoundingBox(3)*imbox(i).BoundingBox(4);
end   
[~, maxsizeind]=max(boxsizes);
imbox = floor(imbox(maxsizeind).BoundingBox);

imbox(imbox<=0) = 1;
width_diff = im_size(2)-(imbox(1)+imbox(3));
if width_diff  < 0 
    imbox(3) = imbox(3)+width_diff;
end
height_diff = im_size(1)-(imbox(2)+imbox(4));
if height_diff  < 0 
    imbox(4) = imbox(4)+height_diff;
end

if any( im_size(1:2) <= roi_size)
    % Our roi size should always be divisible by 2 (for simplicity).
    if rem(max(roi_size),2) ~= 0
        roi_size = max(roi_size)-1;
    end
    
    % If our image is oblong and smaller than a default roi_size, then pad
    % it to be as large as the largest edge.
    padsize = ceil((max(roi_size)-im_size)/2);
    roi = {padarray(test_image,padsize,'both')};
    

else
    % Our roi size should always be divisible by 2 (for simplicity).
    if rem(roi_size,2) ~= 0
        roi_size = roi_size-1;
    end
    roi = cell(round((size(test_image)-roi_size)/roi_step));

    for i=imbox(2):roi_step:imbox(2)+imbox(4)-roi_size
        for j=imbox(1):roi_step:imbox(1)+imbox(3)-roi_size

            numzeros = sum(sum(test_image(i:i+roi_size-1, j:j+roi_size-1)<=10));
            
            if numzeros < (roi_size*roi_size)*0.05
                roi{round(i/roi_step)+1,round(j/roi_step)+1} = test_image(i:i+roi_size-1, j:j+roi_size-1);
            else
                roi{round(i/roi_step)+1,round(j/roi_step)+1} =[];
            end
        end
    end
end

numind = size(roi,1)*size(roi,2);
pixel_spac = nan(size(roi));
confidence = nan(size(roi));

          

% tic;
for r=1:length(pixel_spac(:))
    if ~isempty(roi{r})        
        
        if supersampling% We don't want this run on massive images (RAM sink)

            padsize = roi_size(1)*6; % For reasoning, cite this: https://arxiv.org/pdf/1401.2636.pdf
            padsize = (padsize-roi_size(1))/2 + 1;

            power_spect = fftshift(fft2( padarray(roi{r}, [padsize padsize]) ));
            power_spect = imresize(log10(abs(power_spect).^2),[2048 2048]);
            rhostart = ceil(2048/min(im_size)); % Exclude the DC term from our radial average
        else
            power_spect = fftshift(fft2( roi{r} ));
            power_spect = log10(abs(power_spect).^2);
            rhostart=1; % Exclude the DC term from our radial average
        end
        
%         figure(100); imagesc(roi{r}); axis image;
%         power_spect_export = power_spect-min(power_spect(:));
%         power_spect_export = power_spect_export./max(power_spect_export(:));
%         power_spect_export = power_spect_export.*255;
% 
%         imwrite(uint8(power_spect_export),['pwr_spect ' num2str(r) '.tif']);

        rhosampling = .5;
        thetasampling = 1;

        [polarroi, power_spect_radius] = imcart2pseudopolar(power_spect,rhosampling,thetasampling,[],'makima', rhostart);
        polarroi = circshift(polarroi,-90/thetasampling,1);
        %figure(101); imagesc(polarroi); axis image;
        
        upper_n_lower = [thetasampling:45 136:225 316:360]/thetasampling;
        left_n_right = [46:135 226:315]/thetasampling;
        upper_n_lower_fourierProfile = mean(polarroi(upper_n_lower,:));
        left_n_right_fourierProfile = mean(polarroi(left_n_right,:));
        fullfourierProfile = mean(polarroi);
%         figure(101); plot(upper_n_lower_fourierProfile); axis image;

        if strcmp(row_or_cell,'cell')  && ~all(isinf(left_n_right_fourierProfile)) && ~all(isnan(left_n_right_fourierProfile))

            
            [pixel_spac(r), ~, confidence(r)] = fourierFit(left_n_right_fourierProfile,[], false);
%             [pixel_spac(r), confidence(r)] = fourierFit_rough(left_n_right_fourierProfile, true)
%             confidence(r)
            pixel_spac(r) = 1/ (pixel_spac(r) / ((power_spect_radius*2)/rhosampling));
            
        elseif strcmp(row_or_cell,'row') && ~all(isinf(upper_n_lower_fourierProfile)) && ~all(isnan(upper_n_lower_fourierProfile))

            [pixel_spac(r), ~, confidence(r)] = fourierFit(upper_n_lower_fourierProfile,[], false);
            pixel_spac(r) = 1/ (pixel_spac(r) / ((power_spect_radius*2)/rhosampling));

        else
            pixel_spac(r) = NaN;
        end
    end
    
end
% toc;

avg_pixel_spac = mean(pixel_spac(~isnan(pixel_spac)) );
std_pixel_spac = std(pixel_spac(~isnan(pixel_spac)));
interped_spac_map = avg_pixel_spac;
interped_conf_map = confidence;


%% If we've sampled over the region, then create the heat map
if length(roi) > 1
    interped_spac_map=zeros(im_size);
    interped_conf_map=zeros(im_size);
    interped_corrected_err_map=zeros(im_size);
    sum_map=zeros(im_size);
    
    roi_coverage=roi_size; %round(roi_size/4);
    
    for i=imbox(2):roi_step:imbox(2)+imbox(4)-roi_size
        for j=imbox(1):roi_step:imbox(1)+imbox(3)-roi_size

            if ~isnan( pixel_spac(round(i/roi_step)+1,round(j/roi_step)+1) )
                
                thiserr = confidence(round(i/roi_step)+1,round(j/roi_step)+1)^2;
%                 if thiserr > .44
                    interped_conf_map(i:i+roi_coverage-1, j:j+roi_coverage-1) = interped_conf_map(i:i+roi_coverage-1, j:j+roi_coverage-1) + thiserr;                
                    thisspac = pixel_spac(round(i/roi_step)+1,round(j/roi_step)+1);
                

                    interped_spac_map(i:i+roi_coverage-1, j:j+roi_coverage-1) = interped_spac_map(i:i+roi_coverage-1, j:j+roi_coverage-1) + thiserr*thisspac;

                    sum_map(i:i+roi_coverage-1, j:j+roi_coverage-1) = sum_map(i:i+roi_coverage-1, j:j+roi_coverage-1) + 1;

%                 end
            else

            end
        end
    end
    
    
    
%     [X,Y]=meshgrid( 1:roi_step:(size(test_image,2)-roi_size-1), 1:roi_step:(size(test_image,1)-roi_size-1));
%     [Xq,Yq]=meshgrid( 1:(size(test_image,2)-roi_size-1), 1:(size(test_image,1)-roi_size-1));
%     interped_spac_map = interp2( X,Y, pixel_spac, Xq, Yq);
    
    interped_spac_map = interped_spac_map( imbox(2):imbox(2)+imbox(4), imbox(1):imbox(1)+imbox(3) );
    interped_conf_map = interped_conf_map( imbox(2):imbox(2)+imbox(4), imbox(1):imbox(1)+imbox(3) );
    sum_map = sum_map( imbox(2):imbox(2)+imbox(4), imbox(1):imbox(1)+imbox(3) );
    
    if strcmp(row_or_cell,'cell')
       figure(1);clf; imagesc(interped_spac_map./interped_conf_map); axis image;
    elseif strcmp(row_or_cell,'row')
        figure(1);clf; imagesc((2/sqrt(3)).*interped_spac_map./interped_conf_map); axis image;        
    end
    
    
    
%     [cmap, amap] = firecmap(quantile(scaled_errmap(scaled_errmap~=0), 0.01),...
%                     quantile(scaled_errmap(scaled_errmap~=0), 0.25),...
%                     quantile(scaled_errmap(scaled_errmap~=0), 0.05), 256);
%     [cmap, amap] = firecmap( 0.1882, 0.5608,0.3451, 256);
    
    
    %scaled_confmap = floor(255*(interped_conf_map./sum_map));
    %scaled_confmap(isnan(scaled_confmap)) = 0;
    
    %afullmap = zeros(size(scaled_confmap));
    
    %for i=1:length(afullmap(:))
    %   afullmap(i) = amap( scaled_confmap(i)+1 );
    %end    
    
                
    figure(2);clf; imagesc((interped_conf_map./sum_map)); colormap hot;
%     alpha(errmap,afullmap); axis image;
    
    figure(3);clf; imagesc(sum_map); axis image; colormap gray;
        
%     save( ['Fouriest_Result.mat'], '-v7.3' );
end


% pause;
        
end