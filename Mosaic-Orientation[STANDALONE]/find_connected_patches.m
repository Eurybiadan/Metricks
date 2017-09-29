function [ patches patch_orientations ] = find_connected_patches( im, coords, orientations )
% Robert F Cooper 05-28-2015
% This script finds the connected patches within


triangulation = delaunayTriangulation(coords);
original_orientations = orientations;


alistautocorr=[];
patches = cell(size(coords,1),1);




for i=1:length(triangulation.Points)

    neigh = cell2mat( vertexAttachments(triangulation,i) );
    
    if size(neigh,2) == 6
        
        neighind = unique(triangulation(neigh,:));
        connected_vertices  = coords( neighind,: );
        
        % Determine which vertex consists of the central point
        for v=1:length(connected_vertices)
            if( (connected_vertices(v,1) == triangulation.Points(i,1) ) && ...
                (connected_vertices(v,2) == triangulation.Points(i,2)) )
            
                centerind = neighind(v);

                connected_ind = [neighind(1:v-1,:); neighind(v+1:end,:)];
                break;
            end
        end
        
        % If the central point isn't invalid (-100)
        % Then recursively find connected cells. Otherwise, skip it and
        % move to the next cell in our triangulation
        if orientations(centerind) ~= -100
%             figure(100); 
%             draw_orientation_map( im, coords, original_orientations );
            for ii=1:length(connected_ind)
                
                if orientations(connected_ind(ii)) ~= -100
%                     figure(100); triplot(triangulation);
%                     hold on;                
                    connpatches = vertexAttachments(triangulation,centerind);
%                     triplot( triangulation(connpatches{:},:),triangulation.Points(:,1),triangulation.Points(:,2),'r')
%                     hold off;

                    [conninds, orientations] = conn_patches(triangulation, orientations(centerind), orientations, centerind, connected_ind);


                    % Add the new indexes to the patches list
                    patches{i} = [patches{i}; conninds'];
                                       
                end
            end            
        end
        
    end    
end

% Remove the empty cells
patches = patches(~cellfun('isempty',patches));


% Simplify the regions- if a larger region encompasses a smaller one,
% combine them
for i=1:length(patches)
    
    if ~isempty(patches{i}) && (numel(patches{i}) > 3)
%         i

        % Determine which coordinates are connected to this patch
        patchinds = patches{i};
        boundinds = boundary(coords(patchinds,1), coords(patchinds,2),1);

%         figure(100); 
%         draw_orientation_map( im, coords, original_orientations );
%         plot(coords(patchinds(boundinds),1),coords(patchinds(boundinds),2),'W','LineWidth',5);        
%         hold off;
%         pause(0.5);
        
        % If it can form a boundary, then check to see if it contains any
        % other boundaries or single coordinates.
        if any(boundinds)

            for j=1:length(patches)
                if i ~= j && ~isempty(patches{j})
                    % Determine which vertexes are connected to the potentially internal patch
                    inpatchinds = patches{j};
                    inboundinds = boundary(coords(inpatchinds,1), coords(inpatchinds,2),1);                   
                    
                    if any(inboundinds)
                        [in, on] = inpolygon(coords(inpatchinds(inboundinds),1),coords(inpatchinds(inboundinds),2),...
                                             coords(patchinds(boundinds),1),coords(patchinds(boundinds),2));
                        % Only look at vertexes inside the polygon
                        in = xor(in,on);

                        % If the polygon is contained in the parent polygon,
                        % combine them, and remove the contents of the other
                        % patch.
                        if all(in)
                            figure(100); 
                            draw_orientation_map( im, coords, original_orientations );
                            plot(coords(patchinds(boundinds),1),coords(patchinds(boundinds),2),'g','LineWidth',5);
                            plot(coords(inpatchinds(inboundinds),1),coords(inpatchinds(inboundinds),2),'w','LineWidth',5);
                            hold off;
                            pause(0.5);
                            original_orientations(patches{i})
                            patches{i} = [patches{i}; patches{j}];
                            patches{j} = [];

                            patchinds = patches{i};
                            boundinds = boundary(coords(patchinds,1), coords(patchinds,2),1);
                        end
                    else
                        [in, on] = inpolygon(coords(inpatchinds,1), coords(inpatchinds,2),...
                                             coords(patchinds(boundinds),1),coords(patchinds(boundinds),2));
                                              
                        % Only look at vertexes inside the polygon
                        in = xor(in,on);
                        
                        if all(in)
%                             figure(100); 
%                             draw_orientation_map( im, coords, original_orientations );
%                             plot(coords(patchinds(boundinds),1),coords(patchinds(boundinds),2),'g','LineWidth',5);
%                             plot(coords(inpatchinds,1), coords(inpatchinds,2),'w*','LineWidth',5);
%                             hold off;
%                             pause(0.5);
%                             original_orientations(patches{i})
                            patches{i} = [patches{i}; patches{j}];
                            patches{j} = [];

                            patchinds = patches{i};
                            boundinds = boundary(coords(patchinds,1), coords(patchinds,2),1);
                        end
                        
                    end
                    
                    
                end
            end
        end
    end
end

% Remove the empty patches again
patches = patches(~cellfun('isempty',patches));

patch_orientations = cell(size(patches));

for i=1:length(patches)
   patch_orientations{i} = original_orientations(patches{i});    
end
