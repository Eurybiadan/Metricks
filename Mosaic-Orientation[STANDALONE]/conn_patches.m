function [ conninds, orientations ] = conn_patches( triangulation, ref_orient, orientations, centerind, connected_inds)
% Robert F Cooper 05-28-2015
%   This function looks at the neighboring cells and checks if they are
%   within a tolerance; if they are, they're added to the connected index
%   list and the function is recursively continued.

    TOLERANCE = 3;

if orientations(centerind) ~= -100

    % In the very least, the center index is connected.
    conninds = centerind;

    % Record the orientation
    centero = orientations(centerind);

    % Invalidate the central orientation
    orientations(centerind) = -100;
    tmp_orientations = orientations;
    
    connected_inds= connected_inds(tmp_orientations(connected_inds) ~= -100);   

    if ~isempty(connected_inds)
%         cont = size(tmp_orientations(connected_inds),2)
%         cent = size(centero,2)
        alldists = pdist2(tmp_orientations(connected_inds), centero);

        if any(alldists > 30)

            flipsign = sign(centero - tmp_orientations(connected_inds( alldists  > 30 )));

            tmp_orientations(connected_inds(  alldists > 30 ))  = tmp_orientations(connected_inds( alldists  > 30 )) + flipsign.*60;

    %         tmp_orientations(connected_inds( ( alldists ) <= 30 )) = tmp_orientations(connected_inds( ( alldists ) <= 30 )) + alldists(abs( alldists ) <= 30).*2;
        end



    %     figure(100); hold on;
    %     plot(triangulation.Points(centerind,1),triangulation.Points(centerind,2),'w*','LineWidth',5);
    %     hold off;
    %     pause(0.01);

        for i=1:length(connected_inds)   
            % If the orientation at a neighboring index is valid, then continue
        %     if 
            if tmp_orientations(connected_inds(i)) ~= -100
                % Adjust any phase-wrapping cells to be next to one another

                % Check if the center orientation is within the tolerance
                if abs(tmp_orientations(connected_inds(i))-ref_orient) <= TOLERANCE

                    % Find the neighbors of the cell with the orientation that is
                    % within the tolerance
                    neigh = cell2mat( vertexAttachments(triangulation,connected_inds(i)) );

                    neighind = unique(triangulation(neigh,:));
                    connected_vertices  = triangulation.Points( neighind,: );

                    % Determine which vertex consists of the central point
                    for v=1:length(connected_vertices)
                        if( (connected_vertices(v,1) == triangulation.Points(connected_inds(i),1) ) && ...
                            (connected_vertices(v,2) == triangulation.Points(connected_inds(i),2)) )

                            centerind = neighind(v);
%                             ref_orient = orientations(centerind);
    %                         centerind
                            connected_ind = [neighind(1:v-1,:); neighind(v+1:end,:)];
                            break;
                        end
                    end

    %                 if (ref_orient > 0) && (ref_orient <= 20)
    %                     ref_orient
    %                     orientations(connected_ind)
    %                     goodvert = abs(orientations(connected_ind)-ref_orient) <= TOLERANCE;
    % 
    %                     figure(100); hold on;
    % %                     plot(triangulation.Points(connected_ind(goodvert),1),triangulation.Points(connected_ind(goodvert),2),'g*','LineWidth',5);
    %                     plot(triangulation.Points(centerind,1),triangulation.Points(centerind,2),'w*','LineWidth',5);
    % %                     plot(triangulation.Points(connected_ind(~goodvert),1),triangulation.Points(connected_ind(~goodvert),2),'b*','LineWidth',5);
    %                     hold off;

    %                 end

                    % Recursively call this function
                    % Update the orientation list with any changes we made
                    [conn, orientations] = conn_patches(triangulation,ref_orient, orientations, centerind, connected_ind);
    %                 centerind
                    conninds = [conninds conn];
                else
    %                 ref_orient
    %                 orientations(connected_inds(i))
                end    
            end
        end
    else
        conninds=[];
    end
else
    conninds=[];
end

end

