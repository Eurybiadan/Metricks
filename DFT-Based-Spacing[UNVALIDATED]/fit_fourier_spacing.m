function [avg_pixel_spac, interped_spac_map, interped_err_map, sum_map, imbox ] = fit_fourier_spacing(test_image, roi_size)


if ~exist('test_image','var')
    [filename, pathname] = uigetfile('*.tif', 'Pick an image to segment');

    test_image = imread( fullfile(pathname, filename) );
    if size(test_image,3) >1
        test_image = test_image(:,:,1);
    end
end
% tic;

im_size = size(test_image);
if ~exist('roi_size','var')
    roi_size = 128; %128;%300;
end
roi_step = floor(roi_size/4);
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

if any( im_size < roi_size)    
    roi = {test_image};
else
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
err = nan(size(roi));

  
for r=1:length(pixel_spac(:))
    if ~isempty(roi{r})
        power_spect = fftshift(fft2(roi{r}));
%         all_spect = cat(3,all_spect,abs(power_spect));
        power_spect = log10(abs(power_spect).^2);

        
%         figure(100); imagesc(power_spect); colormap gray; axis image;
%         power_spect_export = power_spect-min(power_spect(:));
%         power_spect_export = power_spect_export./max(power_spect_export(:));
%         power_spect_export = power_spect_export.*255;
% %         
%         imwrite(uint8(power_spect_export),['pwr_spect ' num2str(r) '.tif']);

        rhosampling = .5;
        thetasampling = 1;

        polarroi = imcart2pseudopolar(power_spect,rhosampling,thetasampling,[],'linear');
        polarroi = circshift(polarroi,-90,1);
        
        
        fourierProfile = mean(polarroi);

        if ~all(isinf(fourierProfile)) && ~all(isnan(fourierProfile))

            [pixel_spac(r), ~, err(r)] = fourierFit(fourierProfile,[]);
            pixel_spac(r) = 1/ (pixel_spac(r) / (size(polarroi,2)*2));
            
%             if pixel_spac(r) > 18
%                 fourierFit(fourierProfile,[],true);
%                 smangleProfile = mean(polarroi([1:45 135:225 315:360],:));
%                 figure(1); hold on; plot(smangleProfile-min(smangleProfile)); hold off;
%                 figure(100); imagesc(roi{r}); colormap gray; axis image;
%                 err(r)
%             end
            
        else
            pixel_spac(r) = NaN;
        end

    end
end

avg_pixel_spac = mean(pixel_spac(~isnan(pixel_spac)) );
std_pixel_spac = std(pixel_spac(~isnan(pixel_spac)));
interped_spac_map = avg_pixel_spac;
interped_err_map = err;


%% If we've sampled over the region, then create the heat map
if length(roi) > 1
    interped_spac_map=zeros(im_size);
    interped_err_map=zeros(im_size);
    interped_corrected_err_map=zeros(im_size);
    sum_map=zeros(im_size);
    
    for i=imbox(2):roi_step:imbox(2)+imbox(4)-roi_size
        for j=imbox(1):roi_step:imbox(1)+imbox(3)-roi_size

            if ~isnan( pixel_spac(round(i/roi_step)+1,round(j/roi_step)+1) )
                
                thiserr = err(round(i/roi_step)+1,round(j/roi_step)+1)^2;
%                 if thiserr > .44
                    interped_err_map(i:i+roi_size-1, j:j+roi_size-1) = interped_err_map(i:i+roi_size-1, j:j+roi_size-1) + thiserr;                
                    thisspac = pixel_spac(round(i/roi_step)+1,round(j/roi_step)+1);
                

                    interped_spac_map(i:i+roi_size-1, j:j+roi_size-1) = interped_spac_map(i:i+roi_size-1, j:j+roi_size-1) + thiserr*thisspac;

                    sum_map(i:i+roi_size-1, j:j+roi_size-1) = sum_map(i:i+roi_size-1, j:j+roi_size-1) + 1;

%                 end
            else

            end
        end
    end
    
    
    
%     [X,Y]=meshgrid( 1:roi_step:(size(test_image,2)-roi_size-1), 1:roi_step:(size(test_image,1)-roi_size-1));
%     [Xq,Yq]=meshgrid( 1:(size(test_image,2)-roi_size-1), 1:(size(test_image,1)-roi_size-1));
%     interped_spac_map = interp2( X,Y, pixel_spac, Xq, Yq);
    
    interped_spac_map = interped_spac_map( imbox(2):imbox(2)+imbox(4), imbox(1):imbox(1)+imbox(3) );
    interped_err_map = interped_err_map( imbox(2):imbox(2)+imbox(4), imbox(1):imbox(1)+imbox(3) );
    sum_map = sum_map( imbox(2):imbox(2)+imbox(4), imbox(1):imbox(1)+imbox(3) );
    
    figure(1);clf; imagesc(interped_spac_map./interped_err_map); axis image;
    figure(2);clf; imagesc(interped_err_map./sum_map); axis image; colormap(flipud(jet(256)));
    figure(3);clf; imagesc(sum_map); axis image; colormap gray;
end


% pause;
        
end