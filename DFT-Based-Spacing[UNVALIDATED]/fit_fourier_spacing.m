function [pixel_spac, interped_spac_map] = fit_fourier_spacing(test_image)


if ~exist('test_image','var')
    [filename, pathname] = uigetfile('*.tif', 'Pick an image to segment');

    test_image = imread( fullfile(pathname, filename) );
end
% tic;

im_size = size(test_image);
roi_size = 200;

if any( im_size < roi_size)    
    roi = {test_image};
else
    roi = cell(round(size(test_image)/10)-roi_size);

    for i=1:10:size(test_image,1)-roi_size
        for j=1:10:size(test_image,2)-roi_size

            roi{round(i/10)+1,round(j/10)+1} = test_image(i:i+roi_size-1, j:j+roi_size-1);

        end
    end
end

numind = size(roi,1)*size(roi,2);
pixel_spac = nan(size(roi));

  
for r=1:length(pixel_spac(:))
    if ~isempty(roi{r})
        power_spect = fftshift(fft2(roi{r}));
        power_spect = log10(abs(power_spect).^2);

        rhosampling = .5;
        thetasampling = 1;

        polarroi = imcart2pseudopolar(power_spect,rhosampling,thetasampling,'linear');
        polarroi = circshift(polarroi,-90,1);

        % for i=1:size(polarroi,1)
        %     polarroi(i,:) = conv(polarroi(i,:),[1/5 1/5 1/5 1/5 1/5], 'same' );
        % end

        % figure(1);
        % imagesc(polarroi); colormap gray;

        img = polarroi';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % figure(2);
        % imagesc(power_spect); colormap gray; %hold on; plot(x+size(power_spect,2)/2, y+size(power_spect,1)/2); hold off;

        fourierProfile = mean(polarroi);


        pixel_spac(r) = 1/(fourierFit(fourierProfile,[])); %% THIS IS WRONG!!!
        if pixel_spac(r) > 15
            pixel_spac(r)
        end
%         [i, j]=ind2sub(size(roi),r)
%         spac(r)
    end
end

avg_pixel_spac = mean(pixel_spac);

interped_spac_map=[];
% If we've sampled over the region, then display the heat_map
if length(roi) > 1
    [X,Y]=meshgrid( 1:10:(size(test_image,2)-roi_size-1), 1:10:(size(test_image,1)-roi_size-1));
    [Xq,Yq]=meshgrid( 1:(size(test_image,2)-roi_size-1), 1:(size(test_image,1)-roi_size-1));
    interped_spac_map = interp2( X,Y, pixel_spac, Xq, Yq);
    
%     imagesc(interped_spac_map); axis image;
end



        
end
%% Old, full image approach
% power_spect = fftshift(fft2(test_image));
% power_spect = log10(abs(power_spect).^2);
% 
% rhosampling = .5;
% thetasampling = 1;
% 
% polarroi = imcart2pseudopolar(power_spect,rhosampling,thetasampling,'linear');
% polarroi = circshift(polarroi,-90,1);
% 
% % for i=1:size(polarroi,1)
% %     polarroi(i,:) = conv(polarroi(i,:),[1/5 1/5 1/5 1/5 1/5], 'same' );
% % end
% 
% % figure(1);
% % imagesc(polarroi); colormap gray;
% 
% img = polarroi';
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % figure(2);
% % imagesc(power_spect); colormap gray; %hold on; plot(x+size(power_spect,2)/2, y+size(power_spect,1)/2); hold off;
% 
% fourierProfile = mean(polarroi);
% 
% fourierFit(fourierProfile,[]);
