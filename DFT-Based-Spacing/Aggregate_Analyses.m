clear
close all;

% Find the global dataset
[agg_fName, globalpath ]=uigetfile(fullfile(pwd,'*.mat'),'Select the aggregate file.');

% Individuals will be above this.
individual_path = getparent(globalpath,2,'full');
[fNames] = read_folder_contents(individual_path, 'mat');

% Load the global dataset.
load(fullfile(globalpath, agg_fName));

%% Downsample everything first- we have a lot of redundant data.
downsample_factor = 0.25;

original_size = size(avg_density);
downsampled_size = round(original_size*downsample_factor)+1;
scaling = scaling/downsample_factor;
global_fovea_coords = round(global_fovea_coords*downsample_factor);
montage_rect = cellfun(@(c) round(c*downsample_factor), montage_rect, 'UniformOutput', false);

avg_density = imresize(avg_density, downsampled_size, 'lanczos3');
avg_error = imresize(avg_error, downsampled_size, 'lanczos3');
threshold_mask = imresize(threshold_mask, downsampled_size, 'nearest');
value_std_dev = imresize(value_std_dev, downsampled_size, 'lanczos3');

maskeddensity = avg_density.*threshold_mask;
maskedconf = avg_error.*threshold_mask;


%% Obtain global strip averages of density.
maxstriplen = min([downsampled_size-global_fovea_coords global_fovea_coords]); %Find out the largest strip we can make in pixels.

deg_position = (0:maxstriplen)*scaling;
strip_radius = round(1/(2*scaling)); % Corresponds to 1 degree.
strip_length = length(deg_position)-1;

figure(10); clf; 
% Superior
avg_sup_strip = mean(flipud(maskeddensity(global_fovea_coords(2)-strip_length:global_fovea_coords(2),global_fovea_coords(1)-strip_radius:global_fovea_coords(1)+strip_radius)),2, 'omitnan');
plot(deg_position,avg_sup_strip); hold on;
% Inferior
avg_inf_strip = mean(maskeddensity(global_fovea_coords(2):global_fovea_coords(2)+strip_length, global_fovea_coords(1)-strip_radius:global_fovea_coords(1)+strip_radius),2, 'omitnan');
plot(deg_position,avg_inf_strip);
% Nasal
avg_nasal_strip = mean(maskeddensity(global_fovea_coords(2)-strip_radius:global_fovea_coords(2)+strip_radius,global_fovea_coords(1):global_fovea_coords(1)+strip_length), 1, 'omitnan');
plot(deg_position, avg_nasal_strip);
% Temporal
avg_temp_strip = mean(fliplr(maskeddensity(global_fovea_coords(2)-strip_radius:global_fovea_coords(2)+strip_radius, global_fovea_coords(1)-strip_length:global_fovea_coords(1))), 1, 'omitnan');
plot(deg_position,avg_temp_strip);
legend('Superior','Inferior','Nasal','Temporal')
xlabel('Radial distance (degrees)')
ylabel('Density (cells/degrees^2)')
hold off;

saveas(gcf, fullfile(globalpath,[num2str(length(fNames)) 'subjects_directional_avgdens.svg']) );
saveas(gcf, fullfile(globalpath,[num2str(length(fNames)) 'subjects_directional_avgdens.png']) );

%% Obtain global strip averages of confidence.

figure(11); clf; 
% Superior
avg_sup_strip_conf = mean(flipud(maskedconf(global_fovea_coords(2)-strip_length:global_fovea_coords(2),global_fovea_coords(1)-strip_radius:global_fovea_coords(1)+strip_radius)),2, 'omitnan');
plot(deg_position,avg_sup_strip_conf); hold on;
% Inferior
avg_inf_strip_conf = mean(maskedconf(global_fovea_coords(2):global_fovea_coords(2)+strip_length, global_fovea_coords(1)-strip_radius:global_fovea_coords(1)+strip_radius),2, 'omitnan');
plot(deg_position,avg_inf_strip_conf);
% Nasal
avg_nasal_strip_conf = mean(maskedconf(global_fovea_coords(2)-strip_radius:global_fovea_coords(2)+strip_radius,global_fovea_coords(1):global_fovea_coords(1)+strip_length), 1, 'omitnan');
plot(deg_position, avg_nasal_strip_conf);
% Temporal
avg_temp_strip_conf = mean(fliplr(maskedconf(global_fovea_coords(2)-strip_radius:global_fovea_coords(2)+strip_radius, global_fovea_coords(1)-strip_length:global_fovea_coords(1))), 1, 'omitnan');
plot(deg_position,avg_temp_strip_conf);
legend('Superior','Inferior','Nasal','Temporal')
xlabel('Radial distance (degrees)')
ylabel('Average confidence (AU)')
hold off;

saveas(gcf, fullfile(globalpath,[num2str(length(fNames)) 'subjects_directional_avgconf.svg']) );
saveas(gcf, fullfile(globalpath,[num2str(length(fNames)) 'subjects_directional_avgconf.png']) );

%% Perform a pairwise analysis of the difference of each person from the average strips
figure(12); clf; hold on;
xlabel('Radial distance (degrees)')
ylabel('Density (cells/degrees^2)')
figure(13); clf; hold on;
xlabel('Radial distance (degrees)')
ylabel('Confidence (AU)')

for f=1:length(fNames)
    disp(['Loading, downsampling and prepping: ' fNames{f}])
    load( fullfile(individual_path, fNames{f}), 'density_map_comb', 'blendederr_comb');
    
    indiv_size = [montage_rect{f}(3,2)-montage_rect{f}(1,2)+1 montage_rect{f}(3,1)-montage_rect{f}(1,1)+1];
    indiv_shifted_density = nan(downsampled_size);
    indiv_shifted_confidence = nan(downsampled_size);

    rowrange = round((montage_rect{f}(1,2):montage_rect{f}(3,2))+global_fovea_coords(2));
    colrange = round((montage_rect{f}(1,1):montage_rect{f}(3,1))+global_fovea_coords(1));
    
    if isempty(strfind(fNames{f}, global_eye)) % If it doesn't match our eye, flip the montage data.
        density_map_comb = fliplr(density_map_comb);
        blendederr_comb = fliplr(blendederr_comb);
    end

    
    % Downsample the dataset at the same rate the average was, and mask it.
    density_map_comb = imresize(density_map_comb, indiv_size, 'lanczos3');
    blendederr_comb = imresize(blendederr_comb, indiv_size, 'lanczos3');
   
    
    % Put the dataset it the position we had it before.
    indiv_shifted_density( rowrange, colrange) =  density_map_comb;
    indiv_shifted_confidence( rowrange, colrange) = blendederr_comb;

    indiv_shifted_density = indiv_shifted_density.*threshold_mask;
    indiv_shifted_confidence = indiv_shifted_confidence.*threshold_mask;
    

    indiv_temp_strip = mean(fliplr(indiv_shifted_density(global_fovea_coords(2)-strip_radius:global_fovea_coords(2)+strip_radius,...
                                                         global_fovea_coords(1)-strip_length:global_fovea_coords(1))), 1, 'omitnan');
    figure(12);
    plot(deg_position, avg_temp_strip-indiv_temp_strip);
    drawnow;

    figure(13);
    indiv_temp_strip_conf = mean(fliplr(indiv_shifted_confidence(global_fovea_coords(2)-strip_radius:global_fovea_coords(2)+strip_radius,...
                                                         global_fovea_coords(1)-strip_length:global_fovea_coords(1))), 1, 'omitnan');
    plot(deg_position, avg_temp_strip_conf-indiv_temp_strip_conf);
    drawnow;
end
figure(12); hold off;
figure(13); hold off;