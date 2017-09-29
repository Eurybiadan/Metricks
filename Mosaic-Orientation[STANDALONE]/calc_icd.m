function [ mean_correct_inter_cell_dist ] = calc_icd( coords, um_per_pix )


dt = DelaunayTri(coords);
m=1;
% max_cell_dist=[];
% inter_cell_dist=[];

correct_inter_cell_dist = zeros(size(coords,1),1);
correct_max_cell_dist = zeros(size(coords,1),1);
correct_nn_cell_dist = zeros(size(coords,1),1);
% Find all instances of each coordinate point
for k=1 : size(coords,1)   

    [i j] =find(dt.Triangulation == k);

    conn_ind = dt.Triangulation(i,:);

    coord_row = unique(conn_ind( conn_ind ~= k)); % Find all of the unique coordinate points that isn't the "center" coordinate

    if(size(i,1)~=1)
        coord_row = [k; coord_row]; % Add the "center" to the top, so we know the order for the distances
    else
        coord_row = [k; coord_row']; 
    end

    cell_dist = squareform(pdist([coords(coord_row,1) coords(coord_row,2)]));
            
    % Only take the first row because that is the cell of interest's
    % relative distance to its neighboring cells
    correct_inter_cell_dist(m) = um_per_pix*(sum(cell_dist(1,:)) / (length(cell_dist(1,:))-1));
    correct_max_cell_dist(m)   = um_per_pix*max(cell_dist(1,:));
    correct_nn_cell_dist(m)    = um_per_pix*min(cell_dist(1,2:end));
    m=m+1;

%     inter_cell_dist = [inter_cell_dist scale*(sum(cell_dist(1,:)) / (length(cell_dist(1,:))-1))];
%     max_cell_dist = [max_cell_dist scale*max(cell_dist(1,:))];
end
m = m-1;

% mean_inter_cell_dist = mean(inter_cell_dist);
% mean_max_cell_dist = mean(max_cell_dist);

%     mean_correct_nn_dist = mean( correct_nn_cell_dist(1:m) );
    mean_correct_inter_cell_dist = mean(correct_inter_cell_dist(1:m));
    
%     mean_correct_max_cell_dist   = mean( correct_max_cell_dist(1:m) );

end

