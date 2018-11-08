% Robert F Cooper 2018-11-07
%
%
clear;
close all force;
[fNames,thispath ]=uigetfile(fullfile(pwd,'*.mat'),'Select all montages you wish to combine.', 'MultiSelect', 'on');

foveal_coords = nan(length(fNames), 2);
montage_size = nan(length(fNames), 2);
montage_rect = cell(length(fNames), 1);

for f=1:length(fNames)
    load( fullfile(thispath, fNames{f}), 'fovea_coords', 'imsize');
    
    foveal_coords(f,:) = fovea_coords;
    montage_size(f,:) = imsize(1:2);
    
    montage_rect{f} = [               0                  0; % TLC
                       montage_size(f,2)                 0; % TRC
                       montage_size(f,2) montage_size(f,1); % BRC
                                      0  montage_size(f,1); % BLC
                                      0                  0];% TLC 
    
    % Shift all montage rectangles to a null point
    montage_rect{f} = montage_rect{f}-foveal_coords(f,:);    
end

% Find the bounding rectangle
minglobalbounds = min(cell2mat(montage_rect));
maxglobalbounds = max(cell2mat(montage_rect));

global_dimension = fliplr(ceil(maxglobalbounds-minglobalbounds)); % Flipped so height/width (row/col), not width/height (x/y)

%% Add all of the dft information to the montage, and combine.
combined_spacing = zeros(global_dimension);
combined_error = zeros(global_dimension);
sum_map = zeros(global_dimension);

figure;
for f=1:length(fNames)
    load( fullfile(thispath, fNames{f}), 'fovea_coords', 'imbox', 'im_spac_map',...
                                         'im_err_map', 'im_sum_map');
    
    for i=1:length(imbox)
   
        thisbox = imbox{i}; % x y width height
        thisbox(1:2) = round(thisbox(1:2)-fovea_coords-minglobalbounds);

%         tmprect = [thisbox(1)            thisbox(2); % TLC
%                    thisbox(1)+thisbox(3) thisbox(2); % TRC
%                    thisbox(1)+thisbox(3) thisbox(2)+thisbox(4); % BRC
%                    thisbox(1)            thisbox(2)+thisbox(4); % BLC
%                    thisbox(1)            thisbox(2)]; % TLC
%         plot(tmprect(:,1), tmprect(:,2)); hold on;    

        thismap = im_spac_map{i};
        thiserrmap = im_err_map{i};
        thissummap = im_sum_map{i};

        thismap(isnan(thismap))=0;
        thiserrmap(isnan(thiserrmap))=0;

        combined_spacing( thisbox(2):thisbox(2)+thisbox(4),...
                   thisbox(1):thisbox(1)+thisbox(3) ) = combined_spacing( thisbox(2):thisbox(2)+thisbox(4),...
                                                                   thisbox(1):thisbox(1)+thisbox(3) ) + thismap;

        sum_map( thisbox(2):thisbox(2)+thisbox(4),...
                 thisbox(1):thisbox(1)+thisbox(3) ) = sum_map( thisbox(2):thisbox(2)+thisbox(4),...
                                                               thisbox(1):thisbox(1)+thisbox(3) ) + thissummap;

        combined_error( thisbox(2):thisbox(2)+thisbox(4),...
                   thisbox(1):thisbox(1)+thisbox(3) ) = combined_error( thisbox(2):thisbox(2)+thisbox(4),...
                                                                   thisbox(1):thisbox(1)+thisbox(3) ) + thiserrmap;                                        
    end
end

combined_spacing = combined_spacing./combined_error;
combined_error = combined_error./sum_map;

imagesc(combined_spacing); 
