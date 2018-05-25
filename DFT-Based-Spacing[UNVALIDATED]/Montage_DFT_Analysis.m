clear;
close all;

[thispath]=uigetdir(pwd);

fNames = read_folder_contents(thispath,'tif');

%% Determine the DFT distance for each image in the montage
im_spac_map = cell(length(fNames),1);
im_err_map = cell(length(fNames),1);
im_sum_map = cell(length(fNames),1);
imbox = cell(length(fNames),1);

parfor i=1:length(fNames)
    fNames{i}
    im = imread( fullfile(thispath, fNames{i}) );
    
    if size(im,3) >1
        im = im(:,:,1);
    end
    
    [~, im_spac_map{i}, im_err_map{i}, im_sum_map{i}, imbox{i}] = fit_fourier_spacing(im);    
    
end

%%
imsize = size(imread( fullfile(thispath, fNames{1}) ));

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


blendederrim = blendederrim./sum_map;
blendedim = blendedim./sum_map; % This combination is incorrect because we're averaging averages!


%% Display the results.
figure(1); imagesc(sum_map); colormap gray; axis image;
figure(2); imagesc(blendedim); colormap gray; axis image;
figure(3); imagesc(blendederrim); axis image;

%%
% Empirically determined spacing equation
% to_icd_spac = @(dft_spac) (1.5121*dft_spac.^0.6886);
% 
% scaled_blendedim=to_icd_spac(blendedim.*0.45);
% scaled_blendedim = scaled_blendedim-min(scaled_blendedim(:));
% scaled_blendedim = 255.*(scaled_blendedim./ max(scaled_blendedim(:)));
% 
% figure(4); imagesc( to_icd_spac(blendedim.*0.45) ); axis image;

% To density, assuming perfect packing
row_spac = to_icd_spac(blendedim.*0.45).* (sqrt(3)/ (2));

density_map = (1000^2).*(sqrt(3))./(2 .* row_spac.^2);

figure(5); imagesc( density_map ); axis image;

save( fullfile(thispath,'Fouriest.mat') );
