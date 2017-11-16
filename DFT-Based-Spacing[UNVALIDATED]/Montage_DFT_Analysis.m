clear;
close all;



[thispath]=uigetdir(pwd);

fNames = read_folder_contents(thispath,'tif');

%% Determine the DFT distance for each image in the montage
interped_spac_map = cell(length(fNames),1);
imbox = cell(length(fNames),1);

parfor i=1:length(fNames)
    fNames{i}
    im = imread( fullfile(thispath, fNames{i}) );
    
    if size(im,3) >1
        im = im(:,:,1);
    end
    
    [avg_pixel_spac, interped_spac_map{i}, imbox{i}] = fit_fourier_spacing(im);    
    
end

%%
imsize = size(imread( fullfile(thispath, fNames{1}) ));

blendedim = zeros(imsize);
sum_map = zeros(imsize);

for i=1:length(imbox)
   
    thisbox = imbox{i};
    thismap = interped_spac_map{i};
    thismap(isnan(thismap))=0;
    
    blendedim( thisbox(2):thisbox(2)+thisbox(4),...
               thisbox(1):thisbox(1)+thisbox(3) ) = blendedim( thisbox(2):thisbox(2)+thisbox(4),...
                                                               thisbox(1):thisbox(1)+thisbox(3) ) + thismap;
    
    sum_map( thisbox(2):thisbox(2)+thisbox(4),...
             thisbox(1):thisbox(1)+thisbox(3) ) = sum_map( thisbox(2):thisbox(2)+thisbox(4),...
                                                           thisbox(1):thisbox(1)+thisbox(3) ) + (interped_spac_map{i}>0);
    
end

blendedim = blendedim./sum_map;

figure(1); imagesc(sum_map); colormap gray; axis image;
figure(2); imagesc(blendedim); colormap gray; axis image;