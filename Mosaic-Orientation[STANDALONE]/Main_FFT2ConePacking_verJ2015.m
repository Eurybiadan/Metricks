%==========================================================================
% Giuseppe Lombardo: 24-12-2014
% Requirements: Statistics toolboxes, Image Processing.
%==========================================================================
%% Initialization workspace
% clear
% close all
% clc
%% Add the local path to the search path and select the directory in which 
% are saved.

% mainroutine = 'main_FFT2ConePacking_VerJ2015.m';
% thisPath=which(mainroutine);    % this path
% basePath=thisPath(1:end-length(mainroutine));         % get the absolute path
% 
% % Add the bin directory to run the remainder of the files
% path(path,fullfile(basePath,'bin'))
%% Main constant definition

rotation = 0;    % constant to be used only for reference mosaic coordinates
% um_per_pix = .41;  % scale length micron per pixel - GET THIS FROM FEEDER
% SCRIPT

%% Selection of the directory/file in which are stored the data to be processed
% Select the directory in which are saved the data

% Directory selection

% directory = uigetdir(fullfile(basePath,'data'),'Select the directory in which are saved the data');  
% files = dir(directory);
% numfiles = length(files)-2;

% File selection  MAT FILE in which the simulated coordinates are saved in coord matrix
% [FileNameData,PathNameData,FilterIndexData] = uigetfile('*.mat','Selected cones coordinates file to process');
% load(strcat(PathNameData,FileNameData));

% % File selection  TXT FILE in which the real cone coordinates are saved
% [FileNameData,PathNameData,FilterIndexData] = uigetfile('*.csv','Selected cones coordinates file to process');
% coord = load(fullfile(PathNameData,FileNameData));
% coords = load(fullfile(PathNameData,FileNameData));
coord = coords;

%% Introducing an arbitrary rotation only to the simulated coordinate matrix
% if rotation
%     [theta,rho] = cart2pol(coord(:,1),coord(:,2));
%     
%     rotazione = +7.23*pi/180;                % angle rotation
%     thetanew = theta + rotazione;
%     
%     [xprimo,yprimo] = pol2cart(thetanew,rho);
%     
%     % x0=10;y0=10;width=x_size+200;height=y_size+200;
%     % figure(1)
%     % % set(gcf,'units','points','position',[x0,y0,width,height]);
%     % plot(xprimo,yprimo,'.','MarkerSize',8);
%     % % axis([1 x_size 1 y_size]);
%     % % set(gca,'YDir','reverse');
%     % title('Hexagonal grid');
%     % axis square
%     
%     coord = [xprimo,yprimo];
% end
% 
%% Translation of x & y coordinate in order to have positive coordinates

ax = min(coord(:,1));
ay = min(coord(:,2));

if ax < 0
    coord(:,1) = coord(:,1) - ax + 1;
end
if ay <0
    coord(:,2) = coord(:,2) - ay + 1;
end

ax = max(coord(:,1));
ay = max(coord(:,2));

c = floor(max(coord));

dim1 = nextpow2(c(1));
dim2 = nextpow2(c(2));


if c(1) == 2^dim1
   dim1 = dim1+1; 
end
if c(2) == 2^dim2
   dim2 = dim2+1; 
end

dim = max([dim1,dim2]);

% if (floor(2^dim1/ax))
%     coord(:,1) = coord(:,1) + floor(((2^dim1) - ax)/2)+1;
% end
% if (floor(2^dim2/ay))
%     coord(:,2) = coord(:,2) + floor(((2^dim2) - ay)/2)+1;
% end

if (floor(2^dim/ax))
    coord(:,1) = coord(:,1) + floor(((2^dim) - ax)/2)+1;
end
if (floor(2^dim/ay))
    coord(:,2) = coord(:,2) + floor(((2^dim) - ay)/2)+1;
end

coord = [round(coord(:,1)),round(coord(:,2))];

% %inversion
coord = [coord(:,2),coord(:,1)];

figure(100)
imagesc(im); colormap gray; axis image; hold on;
plot(coord(:,2),coord(:,1),'b.');
hold off;
% title('Hexagonal grid');
% axis square
% axis ij

%     Coord       = load(filename);                       %Prima colonna coordinate Y, Seconda colonna coordinate X;

%% Filling the Image matrix in order to take into account the position of the cones
% dim = max([dim1,dim2]);

Image       = zeros(2^dim,2^dim);
% Image = im;
[X,Y]       = meshgrid(1:1:2^dim);
ind         = sub2ind(size(Image),coord(:,2),coord(:,1));
Image(ind)  = 1;

% Image visualization 
% figure(2)
% subplot(1,2,1);
% imagesc(Image);
% axis tight ij, box on, grid on, axis square

% subplot(1,2,2);
% spy(Image);
% axis tight, box on, grid on, axis square

%% In the case of reading image file
% Image selection 
% [FileName,PathName] = uigetfile('*.tif','Choose the image');
% nomefile            =[PathName,FileName];
% Image               =double(imread(nomefile,'tif'));
% dim                 = nextpow2(size(Image,1));
% [X,Y]= meshgrid(1:1:2^dim);

%%  Extraction of the sub-circular section for computing the cone spacing and direction
% ray = 2^4;                       % ray of the sub-circular region
micronboxsize = mean_nn_dist.*4.5;

ray = ceil(micronboxsize./um_per_pix);

pack = NaN((ceil(size(Image,1)/ray)+1),(ceil(size(Image,2)/ray)+1));
direction = pack;
hexagonsweight = pack;
numfiles = 0;

SPACE = zeros(size(Image));
% PHI = SPACE;
PHI = cell(size(Image));
HEXAG = SPACE;
DIV = SPACE;

% SpaceOrientbar = waitbar(0,'Beginning spacing/orientation calculation ...');

for ic = 1 : ceil(size(Image,2)/ray)
    for ir = 1 : ceil(size(Image,1)/ray)
        numfiles = numfiles+1;
        
%         waitbar(numfiles/(size(pack,1)*size(pack,2)),SpaceOrientbar);
        
        Imagemask = Image;
        xic    = (ic-1)*ray;
        yir     = (ir-1)*ray;
        mask    = (X-xic).^2+(Y-yir).^2 < ray.^2;
        Imagemask(mask==0) = 0;
        
%         figure
%         spy(Imagemask);
        
        [row,col,val]   = find(Imagemask);
        
        if (isempty(row) || (length(row)<6));
            continue
        end
        figure(110);
        plot(xic,yir,'ro');
%         %         coordmask       = [row,col];
 
        rr = row - min(row)+1;
        cc = col - min(col)+1;
        dimsk = nextpow2(max([max(rr),max(cc)]));
        rr = rr + 2^(dimsk-1);
        cc = cc + 2^(dimsk-1);
        
        Imsk           = zeros(2^(dimsk+1),2^(dimsk+1));
        indmask        = sub2ind(size(Imsk),rr,cc);
        Imsk(indmask)  = 1;
        
         figure(101);
         imagesc(Imsk); colormap gray; axis image;
%          imwrite(Imsk.*255,'Image_mask.tif');
%         spy(Imsk);
        
        [spacing,orientation,relativeweight6th]  = FFT2SpacingOrientationExctraction(Imsk,um_per_pix);
         
        pack(ir,ic)             = spacing;          %cone spacing in each sub-region
        direction(ir,ic)        = orientation;      %cone orientation in each sub-region
        hexagonsweight(ir,ic)   = relativeweight6th;%equivalent noise width in each sub-region

        indcal          = sub2ind(size(SPACE),row,col);
        
        SPACE(indcal)   = SPACE(indcal) + spacing;

%         length(indcal)
        
        for ii=1:length(indcal)
            PHI{indcal(ii)}     = [PHI{indcal(ii)}; orientation];
        end

       -PHI{indcal(ii)}.*180/pi;
       
        
%         PHI(indcal)     = PHI(indcal) +orientation;
        HEXAG(indcal)   = HEXAG(indcal) + relativeweight6th;
        DIV(indcal)     = DIV(indcal) + 1;      
        
%         pause(0.1);
%         close all
          
         
    end
end

% waitbar(1,SpaceOrientbar,'Done!');
% pause(1);
% delete(SpaceOrientbar)
%%
% SPACE(ind)  = SPACE(ind)./DIV(ind);
% PHI(ind)    = PHI(ind)./DIV(ind);
% HEXAG(ind)  = HEXAG(ind)./DIV(ind);

for ii=1:length(ind)
    % If we find that the orientation is -30, shift it to postive so
    % that we get a correct mean.
    
    numneg = sum(PHI{ind(ii)} < 0);
    numpos = sum(PHI{ind(ii)} > 0);
    
    if numneg > numpos
        thirty = (round(PHI{ind(ii)}.*180/pi) == 30);
        PHI{ind(ii)}(thirty) = -PHI{ind(ii)}(thirty);    
    else
        thirty = (round(PHI{ind(ii)}.*180/pi) == -30);
        PHI{ind(ii)}(thirty) = -PHI{ind(ii)}(thirty); 
    end
    
    
%     PHI{ind(ii)}
    PHI_new(ind(ii)) = mean(PHI{ind(ii)});
end

PHI = PHI_new';
% Mean results of spacing and orientation inside the photoreceptor mosaic
% cones_spacing   = [nanmean(nanmean(pack)),nanmean(nanstd(pack))];
% cones_direction = [nanmean(nanmean(direction)),nanmean(nanstd(direction))];
% cones_hexag     = [nanmean(nanmean(hexagonsweight)),nanmean(nanstd(hexagonsweight))];


%% Cone spacing and orientation figures

% figure
% % surfc(X,Y,SPACE,'FaceColor','interp','EdgeColor','interp','FaceLighting','phong');
% % axis equal tight, box on, view(2), colorbar;
% surf(X,Y,SPACE,'Parent',axes,'LineStyle','none','FaceColor','interp');
% axis equal tight, box on, view(2), colorbar;
% 
% figure
% % surfc(X,Y,PHI,'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
% % axis tight, box on, grid on, axis square, view(2), colorbar;
% surf(X,Y,PHI,'Parent',axes,'LineStyle','none','FaceColor','interp');
% axis tight, box on, grid on, axis square, view(2), colorbar;

%% Interpolation procedure and representation of the results within the global image

% SPACE_Iterp     = scatteredInterpolant(X(ind),Y(ind),SPACE(ind),'nearest','none');
% PHI_Iterp       = scatteredInterpolant(X(ind),Y(ind),PHI(ind),'nearest','none');
% HEXAG_Iterp     = scatteredInterpolant(X(ind),Y(ind),HEXAG(ind),'nearest','none');
% 
% SPACE_spatial   = SPACE_Iterp(X,Y);
% PHI_spatial     = PHI_Iterp(X,Y);
% HEXAG_spatial   = HEXAG_Iterp(X,Y);
% 
% mask = ~isnan(PHI_spatial);

% figure
% surf(X,Y,SPACE_spatial.*mask,'Parent',axes,'LineStyle','none','FaceColor','interp');
% colormap(jet), caxis([0 max(SPACE_spatial(:))]); colorbar;
% axis equal tight ij, view(2), box on;
% title('Mean hexagonal cones spacing ','FontWeight','bold','FontSize',16);
% set(gca,'FontWeight','bold','FontSize',18);

% Unfold the output. -30 and 30 are actually a 0 rotation, but by default
% it is registering as -30
% peppedegs = (cell2mat(PHI(ind)));
peppedegs = -PHI(ind).*180/pi;

[V,C] = voronoin(coords,{'QJ'});
figure(20);
imagesc(im); colormap gray; axis image; hold on;
hold on;
for i=1:length(C)
   
    vertices=V(C{i},:);
    numedges=size(V(C{i},1),1);
    
    if (all(C{i}~=1)  && all(vertices(:,1)<max(coords(:,1))) && all(vertices(:,2)<max(coords(:,2))) ... % Second row are the Column limits
                     && all(vertices(:,1)>min(coords(:,1))) && all(vertices(:,2)>min(coords(:,2))))
        
%         if ((peppedegs(i) >25) && (peppedegs(i) <= 35)) % 25 to 30, -25 to -30 degs
%             if ((peppedegs(i) >=30) && (peppedegs(i) <= 35))
%                 peppedegs(i) = peppedegs(i)-60; % Adjust to unwrapped orientation
%             end
%             color = 'b';
%         elseif (peppedegs(i) > 35) && (peppedegs(i) <= 45) % -15 to -25 degs
%             peppedegs(i) = peppedegs(i)-60; % Adjust to unwrapped orientation
%             color = 'c';
%         elseif (peppedegs(i) > 45) && (peppedegs(i) <= 55) % -5 to -15 degs
%             color = 'g';
%             peppedegs(i) = peppedegs(i)-60; % Adjust to unwrapped orientation
%         elseif ( (peppedegs(i) >55) && (peppedegs(i) <=61) ) || (peppedegs(i) <= 5)
%             if ( (peppedegs(i) >55) && (peppedegs(i) <=61) ) % -5 to +5 degs
%                 peppedegs(i) = peppedegs(i)-60; % Adjust to unwrapped orientation
%             end
%             color = 'y';
%         elseif (peppedegs(i) >5) && (peppedegs(i) <= 15) % 5 to 15 degs
%             color = 'm';
%         elseif ((peppedegs(i) >15) && (peppedegs(i) <= 25)) % 15 to 25 degs
%             color = 'r';
%         else
%             color = 'k';
%         end 
        
        if ((peppedegs(i) >=-30) && (peppedegs(i) < -25)) || ((peppedegs(i) >25) && (peppedegs(i) <= 30))
            color = 'b';
        elseif (peppedegs(i) >= -25) && (peppedegs(i) < -15)
            color = 'c';
        elseif (peppedegs(i) >=-15) && (peppedegs(i) < -5)
            color = 'g';
        elseif (peppedegs(i) >=-5) && (peppedegs(i) <= 5)
            color = 'y';
        elseif (peppedegs(i) >5) && (peppedegs(i) <= 15)
            color = 'm';
        elseif ((peppedegs(i) >15) && (peppedegs(i) <= 25)) 
            color = 'r';
        else
            peppedegs(i) = -100;
            color = 'k';
        end 

        if (numedges == 6)                
            patch(V(C{i},1), V(C{i},2), ones(size(V(C{i},1))),'FaceColor',color );
        else
            peppedegs(i) = -100;
            patch(V(C{i},1), V(C{i},2), ones(size(V(C{i},1))),'FaceColor','k' );
        end
    else
        peppedegs(i) = -100;
    end
 
end
axis image;
title('Fourier Orientation Map');
hold off;
figure(21); bar(-30:1:30, histc(peppedegs(peppedegs~=-100), -30:1:30),'histc' );

%% Determine the orientation autocorrelation
triangulation = delaunayTriangulation(coords);
peppeautocorr=[];
for i=1:length(triangulation.Points)

    neigh = cell2mat( vertexAttachments(triangulation,i) );
    
    if size(neigh,2) == 6
        
        neighind = unique(triangulation(neigh,:));
        connected_vertices  =coords( neighind,: );
            
        for v=1:length(connected_vertices)
            if( (connected_vertices(v,1) == triangulation.Points(i,1) ) && ...
                (connected_vertices(v,2) == triangulation.Points(i,2)) )

%                 center = connected_vertices(v,:);
                centerind = neighind(v);
%                 connected_vertices = [connected_vertices(1:v-1,:); connected_vertices(v+1:end,:)];
                connected_ind = [neighind(1:v-1,:); neighind(v+1:end,:)];
                break;
            end
        end
        center_orient =peppedegs(centerind);
        if center_orient ~= -100
            conn_orient = peppedegs(connected_ind);
            conn_orient = conn_orient(conn_orient ~=-100);
            peppeautocorr = [peppeautocorr; mean(center_orient-conn_orient)];
        end
    end    
end
figure(22); bar(-30:1:30, histc(peppeautocorr, -30:1:30)./length(peppeautocorr) ,'histc');


% map = [1 0 0;
%        0 1 0;
%        0 0 1;
%        1 1 0;
%        1 0 1;
%        0 1 1];
% map = [1 0 0;
%        0 1 0;
%        0 0 1;
%        1 1 0;
%        1 0 1;
%        0 1 1];
% caxis([0 2.5]);
% colormap(map); colorbar('Ticks',0:.5:2.5, ...
%                         'TickLabels',{'0-10','10-20','20-30','30-40','40-50','50-60'}) ;

% figure
% surf(X,Y, degs,'Parent',axes,'LineStyle','none','FaceColor','interp');
% colormap(jet),  caxis([-30 30]); colorbar; 
% axis equal tight ij, view(2), box on;
title('Peppe Orientation Map','FontWeight','bold','FontSize',16);
% set(gca,'FontWeight','bold','FontSize',18);

% figure
% surf(X,Y,HEXAG_spatial,'Parent',axes,'LineStyle','none','FaceColor','interp');
% colormap(jet), caxis([0 max(HEXAG_spatial(:))]), colorbar;
% axis equal tight ij, view(2), box on;
% title('Equivalent Noise Width','FontWeight','bold','FontSize',16);
% set(gca,'FontWeight','bold','FontSize',18);


%% Routine for calculation index for comparison:

%==========================================================================
% NND evaluation from Rob counting algorithm:
% dist_between_pts = squareform(pdist(coord));                            % Measure the distance from each set of points to the other
% max_ident = eye(length(dist_between_pts)).*max(dist_between_pts(:));    % Make diagonal not the minimum for any observation
% [minval minind] = min(dist_between_pts+max_ident);                      % Find the minimum distance from one set of obs to another
% 
% mean_nn_dist = mean(minval.*um_per_pix); % Distance in microns
% std_nn_dist =  std(minval.*um_per_pix);
% 
% NND = [mean_nn_dist,std_nn_dist];



