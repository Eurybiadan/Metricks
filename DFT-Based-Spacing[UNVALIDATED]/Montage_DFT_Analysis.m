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
    
%     normerrmap = thiserrmap./thissummap;
%     imagesc(normerrmap); pause;
%     thismap(normerrmap<.48) = 0;
%     thissummap(normerrmap<.48)=0;
%     thiserrmap(normerrmap<.48)=0;    
    
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

threshold_mask = (blendederrim>threshold);



scaled_spacing = (blendedim.*0.4567)-min(blendedim(:).*0.4567);
scaled_spacing = 255.*(scaled_spacing./ max(scaled_spacing(:)) );

scaled_error = 255*blendederrim;

allred = ones(256,3);
% allred(:,1) = (255:-1:0)'/255;
% allred(:,2) = (0:255)'/255;
allred(:,3) = (0:255)'/255;


% To density, assuming perfect packing

spac = sqrt(3)./ (2*(blendedim*0.4567).^2);

density_map = (1000^2).*spac;

% To find foveal center
threshdensitymap=density_map;
threshdensitymap(threshdensitymap<90000)=0;
threshdensitymap(isnan(threshdensitymap))=0;

maxes = imregionalmax(threshdensitymap);

s = regionprops(maxes,'centroid');
maxcoords = zeros(length(s),2);
for i=1:length(s)

    maxcoords(i,:)  = s(i).Centroid;
end

fig=figure(10); clf; hold on;
imagesc(density_map); axis image;
plot(maxcoords(:,1),maxcoords(:,2),'b.'); 

ellipsefit = fit_ellipse(maxcoords(:,1),maxcoords(:,2));
[maxcoordsth, maxcoordsr] = cart2pol(maxcoords(:,1)-ellipsefit.X0_in,maxcoords(:,2)-ellipsefit.Y0_in);
f=fit(maxcoordsth,maxcoordsr,'smoothingspline','SmoothingParam',0.9992513623689557);


maxcoordsth = sort(maxcoordsth);
splinefitr = f(maxcoordsth);
[splinefitx, splinefity]= pol2cart(maxcoordsth, splinefitr);
splinefitx = splinefitx + ellipsefit.X0_in;
splinefity = splinefity + ellipsefit.Y0_in;

centerfovea = [ellipsefit.X0_in ellipsefit.Y0_in];


% rotation matrix to rotate the axes with respect to an angle phi
cos_phi = cos( ellipsefit.phi );
sin_phi = sin( ellipsefit.phi );
R = [ cos_phi sin_phi; -sin_phi cos_phi ];

% the ellipse
theta_r         = linspace(0,2*pi);
ellipse_x_r     = ellipsefit.X0 + ellipsefit.a*cos( theta_r );
ellipse_y_r     = ellipsefit.Y0 + ellipsefit.b*sin( theta_r );
rotated_ellipse = R * [ellipse_x_r;ellipse_y_r];

% draw
plot( splinefitx, splinefity,'g.')
plot( rotated_ellipse(1,:),rotated_ellipse(2,:),'r' );

hold off;

foveamask = ~poly2mask(splinefitx,splinefity,size(threshdensitymap,1),size(threshdensitymap,2));

% polardensity = imcart2pseudopolar(density_map, 1, 1, centerfovea);
% figure(11); imagesc(polardensity);
% 
% figure(12); clf; hold on;
% plot(mean([polardensity(1:2,:); polardensity(358:360,:)],'omitnan')); % Nasal
% plot(mean(polardensity(89:93,:),'omitnan')); % Inferior
% plot(mean(polardensity(179:183,:),'omitnan')); % Temporal
% plot(mean(polardensity(269:273,:),'omitnan')); % Superior


scaled_density = density_map-min(density_map(:));
scaled_density = 255.*(scaled_density./ max(scaled_density(:)) );

%% Display and output
figure(1); imagesc(sum_map); axis image; colorbar;

figure(2); imagesc( (blendedim.*0.4567) ); axis image; colorbar;
imwrite( uint8(scaled_spacing.*threshold_mask.*foveamask), parula(256), '11028_OD_thresh_montage_spacing.tif')

figure(3); imagesc(blendederrim); colormap(flipud(jet(256))); axis image; colorbar;
imwrite( uint8(scaled_error), flipud(jet(256)),'11028_OD_thresh_montage_err.tif')

imwrite( uint8(imclose(threshold_mask,ones(11)).*foveamask.*255), '11028_OD_thresh_montage_mask.tif');

figure(5); imagesc( density_map.*imclose(threshold_mask,ones(11)).*foveamask); axis image; colorbar;
imwrite( uint8(scaled_density.*imclose(threshold_mask,ones(11)).*foveamask),parula(256),'11028_OD_thresh_montage_density.tif')

%%
% save( fullfile(thispath,'Fouriest.mat') );
