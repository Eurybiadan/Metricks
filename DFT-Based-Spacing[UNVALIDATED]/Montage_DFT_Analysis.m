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

myPool = parpool();

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

if threshold < 0.3 % Our absolute cutoff for threshold should be 0.3- that is pretty abysmal.
    threshold = 0.3
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
threshdensitymap=density_map.*roi;
quantile(threshdensitymap(threshdensitymap>0),.95)
threshdensitymap(threshdensitymap<quantile(threshdensitymap(threshdensitymap>0),.95))=0;
threshdensitymap(isnan(threshdensitymap))=0;

maxes = imregionalmax(threshdensitymap.*imclose(threshold_mask,ones(11)));

s = regionprops(maxes,'centroid');
maxcoords = zeros(length(s),2);
for i=1:length(s)
    maxcoords(i,:)  = s(i).Centroid;
end

if size(maxcoords,1) > 3
    figure(10); clf; hold on;
    imagesc(density_map.*threshold_mask); axis image;
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

    foveamask = ~poly2mask(splinefitx,splinefity,size(threshdensitymap,1),size(threshdensitymap,2));
else
    warning('Unable to find enough peaks for a foveal estimation.');
    foveamask = ones(size(threshdensitymap));
end
% polardensity = imcart2pseudopolar(density_map, 1, 1, centerfovea);
% figure(11); imagesc(polardensity);
% 
% figure(12); clf; hold on;
% plot(mean([polardensity(1:2,:); polardensity(358:360,:)],'omitnan')); % Nasal
% plot(mean(polardensity(89:93,:),'omitnan')); % Inferior
% plot(mean(polardensity(179:183,:),'omitnan')); % Temporal
% plot(mean(polardensity(269:273,:),'omitnan')); % Superior

%% Display and output
result_path = fullfile(thispath,'Results');
mkdir(result_path)

fnamesplits = strsplit(fNames{1},'_');
prefix=[];
for f=1:5 % build our prefix.
    prefix = [prefix fnamesplits{f} '_'];
end

%%
save( fullfile(result_path,[prefix 'Fouriest_Result.mat']), '-v7.3' );

%% Show figures

% figure(1); imagesc(sum_map); axis image; colorbar;

blendedim = blendedim.*2/sqrt(3); % To convert to ICD spacing.
scaled_spacing = (blendedim.*scaling)-min(blendedim(:).*scaling);
scaled_spacing = 255.*(scaled_spacing./ max(scaled_spacing(:)) );

figure(2); imagesc( (blendedim.*scaling) ); axis image; colorbar;
imwrite( uint8(scaled_spacing.*threshold_mask.*foveamask), parula(256), fullfile(result_path, [prefix 'thresh_montage_spacing_' num2str(scaling,'%5.2f') '.tif']))
clear scaled_spacing;

scaled_error = 255*blendederrim;
figure(3); imagesc(blendederrim); colormap(flipud(jet(256))); axis image; colorbar;
imwrite( uint8(scaled_error), flipud(jet(256)),fullfile(result_path, [prefix 'thresh_montage_err_' num2str(scaling,'%5.2f') '.tif']))
clear scaled_error;

imwrite( uint8(imclose(threshold_mask,ones(11)).*foveamask.*255), fullfile(result_path, [prefix 'thresh_montage_mask_' num2str(scaling,'%5.2f') '.tif']));

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

