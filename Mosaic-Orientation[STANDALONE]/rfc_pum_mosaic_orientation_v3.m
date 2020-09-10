
%% Robert Cooper 03-07-2014
% This script reimplements the work by Pum et. al in 1990...

% clear;
% close all;

[fname pathname] = uigetfile('*.txt;*.csv');

% coords = dlmread('RFC_68600_cones_mm_0p39486_coords.csv');
% coords = dlmread('RFC_68600_cones_mm_0p39486_coords_jittered.csv');
% coords = dlmread('RFC_68600_cones_mm_0p39486_coords_-30deg.csv');
coords = dlmread(fullfile(pathname,fname));
im = imread(fullfile(pathname, [fname(1:end-11) '.tif'] ));
% coords = dlmread('angle_phantom.csv');

triangulation = delaunayTriangulation(coords);
all_deg=ones(size(coords,1),1).*-100;
all_centers=[];
degs=[];
direction=[];
all_angle=[];
all_sides=[];

for i=1:length(triangulation.Points)
    
%     if bound(i)
        
        neigh = cell2mat( vertexAttachments(triangulation,i) );

        total_angle = (size(neigh,2)-2)*180;
        range = (total_angle/size(neigh,2))/2; % Half of one of the interior angles
        
        all_sides=[all_sides; size(neigh,2)];
%         if size(neigh,2) == 6
            
            connected_vertices  =coords(unique(triangulation(neigh,:)),:);
            
            for v=1:length(connected_vertices)
                if( (connected_vertices(v,1) == triangulation.Points(i,1) ) && ...
                    (connected_vertices(v,2) == triangulation.Points(i,2)) )
                    
                    center = connected_vertices(v,:);
                    connected_vertices = [connected_vertices(1:v-1,:); connected_vertices(v+1:end,:)];
                    break;
                end
            end
            
%             figure(1);
%             triplot(triangulation);
%             hold on;
%             triplot(triangulation(neigh,:),coords(:,1),coords(:,2),'Color','r'); 
%             hold off;
%             drawnow;
            
            % Make vectors out of all of the neighbors
            neighbor_vectors = connected_vertices - repmat(center, size(connected_vertices,1), 1);
            
            % Convert them to magnitudes
            mag_vectors = sqrt( sum(neighbor_vectors.^2, 2) );
            % Find the median- assume the 4th element would be the median
            % (any 7th point would be further away)
            [medianvect medianind] = sort(mag_vectors);
            % Use it as the reference vector (that all will be backrotated
            % to)
            ref_vector  = neighbor_vectors(medianind(4),:);
%             neighbor_vectors = neighbor_vectors(medianind([1:3 5:6]),:);
            
            % Normalize the vectors
            norm_ref_vector  = ref_vector./mag_vectors(medianind(4));            
            ref_angle        = atan2d(norm_ref_vector(2), norm_ref_vector(1));
            
            if( ref_angle < 0 )
               ref_angle = 360+ref_angle; 
            end
            
            while ref_angle > (range/2)
               ref_angle = ref_angle-range;
            end
                         
            norm_backrot_vector = [];
            all_angle           = [];

            norm_neighbor_vectors = neighbor_vectors./repmat(mag_vectors,1,2);

%             if size(neighbor_vectors, 1) == 7
%                 size(neighbor_vectors)
%             end 
            clear back_rot_vector;
            for j=1:size(norm_neighbor_vectors, 1)
                
                angle        = atan2d(norm_neighbor_vectors(j,2), norm_neighbor_vectors(j,1));
            
                if( angle < 0 ) % If it's negative, add 360 so we're in a 0-360 range.
                   angle = 360+angle; 
                end
                
                
%                 while angle > 30
%                    angle = angle-60;
%                 end

                while abs(angle-ref_angle) > (range/2)
                   angle = angle- sign(angle-ref_angle)*range;
                end
                                
                back_rot_vector(j,:) = [ cosd(angle), sind(angle) ];

                all_angle(j) = angle;

            end

            
            tmp = mean( back_rot_vector);
            
            tmpmag = sqrt( sum(tmp.^2) );
            tmp = tmp./tmpmag;
            
%             direction = [direction ; sum(norm_backrot_vector)./6];
            degs      = [degs; atan2d(tmp(2), tmp(1)) ];
            direction = [direction ; tmp];
            all_centers = [all_centers; center];
            
            if abs(degs(end)) > (range/2)
                degs(end) = degs(end) - sign(degs(end)).*range;
            end
            all_deg(i) = -(degs(end));
            
            clear norm_backrot_vector
            clear norm_backrot_vector_offset
%             clear norm_ref_vector
            clear all_angle
           
            
            
%         end
        
%     end
    
end


% triplot(delaunay(coords(:,1), coords(:,2)), coords(:,1), coords(:,2));

% h1 = axes('Parent',gcf,'Color',[0 0 0]);
% hold on;
% voronoi(coords(:,1), coords(:,2), 'k');
% hold on;

%%

% im = imread(fullfile(pathname, [fname(1:end-11) '.tif'] )); 
figure;
imagesc(im); colormap gray; axis image; 
hold on;

% colorgrad = jet( max(six_sided_deg)+1 );

% direction2 = [cosd(tempoutput(:,4)) sind(tempoutput(:,4))];
% all_centers2 = tempoutput(:,1:2);
% ori_region = zeros([size(im) 3]);
% ori_conv = zeros([size(im) 3]);
% 
for i = 1:size(all_centers,1)
    
    switch(all_sides(i))
        case 4
            color = 'm';
        case 5
            color = 'c';
        case 6
            color = 'g';

        case 7
            color = 'y';
        case 8
            color = 'r';
        case 9
            color = 'b';
        otherwise
            color = 'k';
    end



    quiver(all_centers2(i,1), all_centers2(i,2), direction2(i,1), direction2(i,2), 5,'Color',color,'LineWidth',2);
%     plot(six_sided(i,1), six_sided(i,2),'Color',colorgrad(degs(i).*3,:),'Marker','*','MarkerSize',5);
    
end
axis image;
hold off;
return
%%
% colorgrad = jet(6);
% colormap(colorgrad);
% caxis([-30 30])
[V,C] = voronoin(coords,{'QJ'});
figure(10);
% im = imread(fullfile(pathname, [fname(1:end-11) '.tif'] ));
imagesc(im); axis image; colormap gray;hold on;
% distribut = [];
% Color the voronoi cells based on the above values
for i=1:length(C)
   
    vertices=V(C{i},:);
    numedges=size(V(C{i},1),1);
    
    if (all(C{i}~=1)  && all(vertices(:,1)<max(coords(:,1))) && all(vertices(:,2)<max(coords(:,2))) ... % Second row are the Column limits
                      && all(vertices(:,1)>min(coords(:,1))) && all(vertices(:,2)>min(coords(:,2))))
        
        if ((all_deg(i) >=-30) && (all_deg(i) < -25)) || ((all_deg(i) >25) && (all_deg(i) <= 30))
            color = 'b';
        elseif (all_deg(i) >= -25) && (all_deg(i) < -15)
            color = 'c';
        elseif (all_deg(i) >=-15) && (all_deg(i) < -5)
            color = 'g';
        elseif (all_deg(i) >=-5) && (all_deg(i) <= 5)
            color = 'y';
        elseif (all_deg(i) >5) && (all_deg(i) <= 15)
            color = 'm';
        elseif ((all_deg(i) >15) && (all_deg(i) <= 25)) 
            color = 'r';
        else
            color = 'k';
        end 

%         if i == 442
%             i
%         end
        
        if (numedges == 6)                
            patch(V(C{i},1), V(C{i},2), ones(size(V(C{i},1))),'FaceColor',color );
        else
            all_deg(i) = -100;
            patch(V(C{i},1), V(C{i},2), ones(size(V(C{i},1))),'FaceColor','k' );
        end
    else
%         if i == 442
%             i
%         end
        
        all_deg(i) = -100;
    end

end
axis image;
title('Pum Orientation Map')
hold off;
figure(11); bar(-30:1:30, histc(all_deg(all_deg~=-100), -30:1:30),'histc' );
title('Pum Orientation Histogram'); xlabel('Orientation (degrees)'); ylabel('Number of Cells');

%% Determine the orientation autocorrelation
triangulation = delaunayTriangulation(coords);
six_sidedautocorr=[];
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
        center_orient =all_deg(centerind);
        if center_orient ~= -100
            conn_orient = all_deg(connected_ind);
            conn_orient = conn_orient(conn_orient ~=-100);
            six_sidedautocorr = [six_sidedautocorr; mean(center_orient-conn_orient)];
        end
    end    
end
figure(12); bar(-30:1:30, histc(six_sidedautocorr, -30:1:30)./length(six_sidedautocorr) ,'histc');
