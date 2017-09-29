clear;
close all force; 

[fname pathname] = uigetfile('*.txt;*.csv');

coords = dlmread(fullfile(pathname,fname));

coords = coords;

imname = [fname(1:end-11) '.tif'];
impath = fullfile(pathname,imname);
if exist(impath,'file')
   im = imread(impath);
end
   
um_per_pix=.5;

% Calculate the cell distances
mean_nn_dist=calc_icd(coords,um_per_pix);

% figure(33); imagesc(im); colormap gray; axis image;
% recth = rectangle('Position',[0 0 1 1], 'EdgeColor','r');

% micronboxsize = 12;
micronboxsize = mean_nn_dist.*4.5;

pixelboxsize = ceil(micronboxsize./um_per_pix);
pixelboxhalf = floor(pixelboxsize/2);

% center it on 0 so that we can rotate/translate to the center

anglecoords = [-3, 0;
                3, 0];

[X, Y]=meshgrid(0:(pixelboxhalf*2), 0:(pixelboxhalf*2));

X = X-pixelboxhalf;
Y = Y-pixelboxhalf;

mask = uint8(X.^2 + Y.^2 < pixelboxhalf.^2);


[V,C] = voronoin(coords,{'QJ'});

% Determine the number of sides for each coordinate location
for i=1:length(C)
   
    vertices=V(C{i},:);
    
    
if (all(C{i}~=1)  && all( vertices(:,1)<=max(coords(:,1))) && all(vertices(:,2)<=max(coords(:,2)) ) ... % Second row are the Column limits
                  && all( vertices(:,1)>=min(coords(:,1))) && all(vertices(:,2)>=min(coords(:,2)) ) )

    numedges(i)=size(V(C{i},1),1);

end
 
end



i=1;
j=1;
alist=-100*ones(size(coords,1),1);
for i=1:size(coords,1)
        
    x = round(coords(i,1));
    y = round(coords(i,2));
    
    if (y-pixelboxhalf) > 0           && (x-pixelboxhalf) > 0 && ...
       (y+pixelboxhalf) <= size(im,1) && (x+pixelboxhalf) <= size(im,2) ...
       && numedges(i)~=0
       
        roi = im( y-pixelboxhalf : y+pixelboxhalf, ...
                  x-pixelboxhalf : x+pixelboxhalf);

       
        roi = roi.*mask;
        outroi = double(roi.*mask);

        expected_angle(i) = (180*(numedges(i)-2))/numedges(i);
        range{i} = (90 - (expected_angle(i)/4)) : (90 + (expected_angle(i)/4));
        
        radonim = radon(roi,range{i})';

        clear numpeaks;
        
        % Determine the cutoffs for the radon xform by looking for the FWHM
        gausfilt = fspecial('gaussian',[1 5],.75);
        middle_orient = round( size(radonim,1)/2);
        
        middlerow = conv(radonim(middle_orient,:),gausfilt,'valid');
        
        xings = 1:length(middlerow);
        xings = xings(middlerow > (max(radonim(middle_orient,:))/2) );
%         xings = find_zero_crossings(diff(middlerow,2));
         
        rrms=[];
        for r=1:size(radonim,1)
            radonrow = radonim(r,:);
            
            
            radonrow = conv(radonrow,gausfilt,'valid');
            radonrow = radonrow( (xings(1)):(xings(end)));
            dervline =  ( diff( radonrow,2) ); %, 4 );
            
            rrms(r) = rms( dervline );


%             pause(1);
        end 
        
        rrms= rrms';
%         if any(numpeaks>=3) % For edge cases in generated mosaics...
%             numpeaks = numpeaks>=3;
%             rrms = rrms(numpeaks);
%         end
        
                  
        [rmsval angle] = max(rrms);
        [val worstangle] = min(rrms);
        
        theta = (angle - (1+(expected_angle(i)/4)) );

        alist(i) = theta;
        arms(i) = rmsval;
    end
end


    
% figure(1); imagesc(im); axis image; colormap gray;
figure(30); imagesc(im); colormap gray; axis image; hold on;

[V,C] = voronoin(coords,{'QJ'});

green=0;
blue=0;

anglecoords = [-3, 0;
                3, 0];

for ii=1:length(coords)

    x = (coords(ii,1));
    y = (coords(ii,2));

    if alist(ii) ~= -100
        % Draw the orientations of each mosaic on the image
        rfctform = affine2d([cosd(alist(ii)) -sind(alist(ii)) 0; sind(alist(ii)) cosd(alist(ii)) 0; x y 1]);

        [rfcrotx rfcroty]=transformPointsForward(rfctform,anglecoords(:,1),anglecoords(:,2));

        switch(numedges(ii))
            case 4
                linecolor = 'm';
            case 5
                linecolor = 'c';
            case 6
                linecolor = 'g';
            case 7
                linecolor = 'y';
            case 8
                linecolor = 'r';
            case 9
                linecolor = 'b';            
        end
        

        plot(rfcrotx,rfcroty,linecolor); 
        linecolor='w';
        
    end

end
hold off;
figure(30);
% saveas(gcf, [fname(1:end-4) '_radon_allsides.eps'],'epsc');
saveas(gcf, [fname(1:end-4) '_radon_allsides.tif']);

figure(31); imagesc(im); colormap gray; axis image; hold on;
%% Color the voronoi cells based on the above values
for i=1:length(C)
   
    vertices=V(C{i},:);
    cmap = parula(length(range{i}));
    
if (all(C{i}~=1)  && all( vertices(:,1)<=max(coords(:,1))) && all(vertices(:,2)<=max(coords(:,2)) ) ... % Second row are the Column limits
                  && all( vertices(:,1)>=min(coords(:,1))) && all(vertices(:,2)>=min(coords(:,2)) ) ) ...
                  && alist(i) ~= -100
        
 
%         if ((alist(i) >=-30) && (alist(i) < -25))
% %             if alist(i) == -30 
% %                 alist(i) = -100;
% %                 color = 'k';
% %             else
% %              
%                 color = 'b';
% %             end
%         elseif (alist(i) >= -25) && (alist(i) < -15)
%             color = 'c';
%         elseif (alist(i) >=-15) && (alist(i) < -5)
%             color = 'g';
%             green = green+1;
%         elseif (alist(i) >=-5) && (alist(i) <= 5)
%             color = 'y';
%         elseif (alist(i) >5) && (alist(i) <= 15)
%             color = 'm';
%         elseif ((alist(i) >15) && (alist(i) <= 25)) 
%             color = 'r';
%         elseif  ((alist(i) >25) && (alist(i) <= 30))
%             color = 'b';
%         else
% %             alist(i) = -100;
%             color = 'k';
%         end 
        if alist(i) ~= -100
            c = cmap( round(alist(i) + (expected_angle(i)/4)+1 ),:);
        else
            patch(V(C{i},1), V(C{i},2), ones(size(V(C{i},1))),'FaceColor','k' );
        end
        
        if (numedges(i) ~= 0)            
            patch(V(C{i},1), V(C{i},2), ones(size(V(C{i},1))),'FaceColor',c );
%         else
%             alist(i) = -100;
%             patch(V(C{i},1), V(C{i},2), ones(size(V(C{i},1))),'FaceColor','k' );
        end
end
 
end
axis image;
title('Radon Orientation Map')
colorbar;
hold off;

% figure(31); bar(-30.5:1:30.5,histc(alist(alist~=-100), -30.5:1:30.5) ,'histc');


%% Determine the orientation autocorrelation
% triangulation = delaunayTriangulation(coords);
% alistautocorr=[];
% for i=1:length(triangulation.Points)
% 
%     neigh = cell2mat( vertexAttachments(triangulation,i) );
%     
%     if size(neigh,2) == 6
%         
%         neighind = unique(triangulation(neigh,:));
%         connected_vertices  =coords( neighind,: );
%             
%         for v=1:length(connected_vertices)
%             if( (connected_vertices(v,1) == triangulation.Points(i,1) ) && ...
%                 (connected_vertices(v,2) == triangulation.Points(i,2)) )
% 
% %                 center = connected_vertices(v,:);
%                 centerind = neighind(v);
% %                 connected_vertices = [connected_vertices(1:v-1,:); connected_vertices(v+1:end,:)];
%                 connected_ind = [neighind(1:v-1,:); neighind(v+1:end,:)];
%                 break;
%             end
%         end
%         
%         if alist(centerind) ~= -100
%             conn_orient = alist(connected_ind);
%             conn_orient = conn_orient(conn_orient ~=-100);
%             alistautocorr = [alistautocorr; mean(alist(centerind)-conn_orient)];
%         end
%     end    
% end
% figure(32); bar(-30:1:30,histc(alistautocorr, -30:1:30)./length(alistautocorr) ,'histc');