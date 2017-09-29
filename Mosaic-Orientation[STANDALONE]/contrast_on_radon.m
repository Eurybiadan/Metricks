
clear;
im = imread('/remote_project_folders/Rob Cooper/Cell_Orientation_2015/AOSLO_Images/JC_10145_790nm_OD_confocal_0018_ref_44_lps_12_lbss_12_sr_n_70_cropped_5_100um T.tif');
micronboxsize = 12;
mppix = .4;

pixelboxsize = round(micronboxsize./mppix);
pixelboxhalf = round(pixelboxsize/2);


[X, Y]=meshgrid(0:pixelboxsize, 0:pixelboxsize);

X = X-pixelboxhalf;
Y = Y-pixelboxhalf;

mask = uint8(X.^2 + Y.^2 < pixelboxhalf.^2);

offset = ceil(pixelboxsize/2)+1;

i=1;
j=1;
alist=-ones(size(im));
for i=offset:size(im,1)-offset

    for j=offset:size(im,2)-offset
    
    roi = histeq( im( i-pixelboxhalf : i+pixelboxhalf, ...
                      j-pixelboxhalf : j+pixelboxhalf) );

    roi = roi.*mask;
                      
        radonim = radon(roi,60:120)';

%     figure(1); imagesc(roi); colormap gray; axis image;
%     figure(2); imagesc(radonim(:,10:30)); colormap gray; axis image;

    for r=1:size(radonim,1)

        rms(r) = std( diff(radonim(r, 10:30)));

    end 
    rms= rms';

    [val angle]=max(rms);


%     figure(3); plot(radonim(angle,:));
%     figure(4); plot(diff(radonim(angle,:)) );

%     disp(['Found angle was: ' num2str(55+angle)]);
    
    alist(i,j) = 60-angle;
    
%     pause(.5);
    end
end


abinned = zeros(size(alist,1), size(alist,2), 3);

for i=1:size(alist,1)
    for j=1:size(alist,2)
   
        if (alist(i,j) >=0) && (alist(i,j) < 10)
            abinned(i,j,:) = [1 0 0];
        elseif (alist(i,j) >= 10) && (alist(i,j) < 20)
            abinned(i,j,:) = [0 1 0];
        elseif (alist(i,j) >=20) && (alist(i,j) < 30)
            abinned(i,j,:) = [0 0 1];
        elseif (alist(i,j) >=30) && (alist(i,j) < 40)
            abinned(i,j,:) = [1 1 0];
        elseif (alist(i,j) >=40) && (alist(i,j) < 50)
            abinned(i,j,:) = [1 0 1];
        elseif (alist(i,j) >=50) && (alist(i,j) <= 60)
            abinned(i,j,:) = [0 1 1];
        else
            abinned(i,j,:) = [0 0 0];
        end 
        
    end
end

figure(30); imagesc(abinned); axis image;
title('Radon Orientation Map- continuous')