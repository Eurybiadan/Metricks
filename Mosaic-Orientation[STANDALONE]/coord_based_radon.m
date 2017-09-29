
clear;

[fname pathname] = uigetfile('*.txt;*.csv');

coords = dlmread(fullfile(pathname,fname));
im = imread(fullfile(pathname,[fname(1:end-11) '.tif']) );

% im = imread('/remote_project_folders/Rob Cooper/Cell_Orientation_2015/AOSLO_Images/JC_10145_790nm_OD_confocal_0018_ref_44_lps_12_lbss_12_sr_n_70_cropped_5_100um T.tif');
% coords = dlmread('/remote_project_folders/Rob Cooper/Cell_Orientation_2015/AOSLO_Images/JC_10145_790nm_OD_confocal_0018_ref_44_lps_12_lbss_12_sr_n_70_cropped_5_100um T_coords.txt');
micronboxsize = 12;
mppix = .41;

pixelboxsize = round(micronboxsize./mppix);
pixelboxhalf = round(pixelboxsize/2);


[X, Y]=meshgrid(0:pixelboxsize, 0:pixelboxsize);

X = X-pixelboxhalf;
Y = Y-pixelboxhalf;

mask = uint8(X.^2 + Y.^2 < pixelboxhalf.^2);


% figure(31); voronoi( coords(:,1), coords(:,2) ); axis image;

i=1;
j=1;
alist=-ones(size(im));
for i=1:size(coords,1)
        
    x = round(coords(i,1));
    y = round(coords(i,2));
    
    if (y-pixelboxhalf) > 0 && (x-pixelboxhalf) > 0 && ...
       (y+pixelboxhalf) < size(im,1) && (x+pixelboxhalf) < size(im,2)
       
        roi = histeq( im( y-pixelboxhalf : y+pixelboxhalf, ...
                          x-pixelboxhalf : x+pixelboxhalf) );

        roi = roi.*mask;
                      
        radonim = radon(roi,60:120)';

%         radonim = radonim(,:);

%         figure(1); imagesc(roi); colormap gray; axis image;
%         figure(2); imagesc(radonim); colormap gray; axis image;

        for r=1:size(radonim,1)

            dervline = diff(radonim(r,10:38));
            
            rms(r) = std( dervline );
            stripcontrast(r) = (max(radonim(r,10:38)) - min(radonim(r,10:38))) ./ (max(radonim(r,10:38)) + min(radonim(r,10:38)));
            
%             fftim(r,:) = abs(fft(radonim(r,:),[],2) );
            
            
        end 
        rms= rms';

%         figure(30); imagesc(fftim); axis image;
%         pause;
        
        [val angle]=max(rms);
        [val anglecont]=max(stripcontrast);

        disp( [num2str(angle) '/' num2str(anglecont)] );

    %     figure(3); plot(radonim(angle,:));
    %     figure(4); plot(diff(radonim(angle,:)) );

    %     disp(['Found angle was: ' num2str(55+angle)]);

        alist(i) = 60-angle;
%         pause(.4);
    end
end

% figure(1); imagesc(im); axis image; colormap gray;
figure(30);
[V,C] = voronoin(coords,{'QJ'});

% Color the voronoi cells based on the above values
for i=1:length(C)
   
    vertices=V(C{i},:);
    numedges=size(V(C{i},1),1);
    
if (all(C{i}~=1)  && all( vertices(:,1)<max(coords(:,1))) && all(vertices(:,2)<max(coords(:,2)) ) ... % Second row are the Column limits
                  && all( vertices(:,1)>min(coords(:,1))) && all(vertices(:,2)>min(coords(:,2)) ) )
        
        if (alist(i) >=0) && (alist(i) < 10)
            color = 'r';
        elseif (alist(i) >= 10) && (alist(i) < 20)
            color = 'g';
        elseif (alist(i) >=20) && (alist(i) < 30)
            color = 'b';
        elseif (alist(i) >=30) && (alist(i) < 40)
            color = 'y';
        elseif (alist(i) >=40) && (alist(i) < 50)
            color = 'm';
        elseif (alist(i) >=50) && (alist(i) <= 60)
            color = 'c';
        else
            color = 'k';
        end 
        if (numedges == 6)
            patch(V(C{i},1), V(C{i},2), ones(size(V(C{i},1))),'FaceColor',color );
        else
            patch(V(C{i},1), V(C{i},2), ones(size(V(C{i},1))),'FaceColor','k' );
        end
end
 
end
axis image;
% title('Radon Orientation Map')
hold off;

