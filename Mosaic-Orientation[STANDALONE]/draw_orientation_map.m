function [ ] = draw_orientation_map( im, coords, orientations )
%Robert F Cooper 05-28-2015
% This function draws an orientation map given coords, orientations, and an
% image

figure(100); imagesc(im); colormap gray; axis image; hold on;

[V,C] = voronoin(coords,{'QJ'});

for i=1:length(C)
   
    vertices=V(C{i},:);
    numedges=size(V(C{i},1),1);
    
if (all(C{i}~=1)  && all( vertices(:,1)<=max(coords(:,1))) && all(vertices(:,2)<=max(coords(:,2)) ) ... % Second row are the Column limits
                  && all( vertices(:,1)>=min(coords(:,1))) && all(vertices(:,2)>=min(coords(:,2)) ) )
        
 
        if ((orientations(i) >=-30) && (orientations(i) < -25))
            color = 'b';
        elseif (orientations(i) >= -25) && (orientations(i) < -15)
            color = 'c';
        elseif (orientations(i) >=-15) && (orientations(i) < -5)
            color = 'g';            
        elseif (orientations(i) >=-5) && (orientations(i) <= 5)
            color = 'y';
        elseif (orientations(i) >5) && (orientations(i) <= 15)
            color = 'm';
        elseif ((orientations(i) >15) && (orientations(i) <= 25)) 
            color = 'r';
        elseif  ((orientations(i) >25) && (orientations(i) <= 30))
            color = 'b';
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
title('Orientation Map')
% hold off;


end

