%% Robert F Cooper
% 1/13/21
%
% This script's sole job is to run multiple montages in series; to save time, 
% instead of operating on individual files though, it operates on a series
% of folders.

folderList = {};

thisfolder = pwd;



thisfolder = uigetdir(thisfolder, 'Select the folders containing the montages you wish to assess.');
    
[folderList, isdir]=read_folder_contents(thisfolder);

[lutfname, lutfolder] = uigetfile(fullfile(pwd,'*.csv'),'Select scaling LUT, OR cancel if you want to input the scale directly.');

%%
restartf=1;
% endf=1;
endf=length(folderList);
%%
for f=restartf:endf
    restartf=f;
    confocal_coords=[];
    something=false;
    if isdir(f)
        confocaldir=fullfile(thisfolder, folderList{f},'confocal');
        splitdir=fullfile(thisfolder, folderList{f},'split detection');
        
        if exist(confocaldir, 'dir')
            disp(['***** Running analysis on: ' confocaldir ' *****']);
            fNames = read_folder_contents(confocaldir,'tif');

            [ scalinginfo, ~, lut ]=determine_scaling(confocaldir, fNames, fullfile(lutfolder,lutfname) ,'degrees');

            confocal_coords=Foveated_Montage_DFT_Analysis(confocaldir, fNames, scalinginfo, 'degrees', lut, true);
            something=true;
        end
        
        if exist(splitdir, 'dir')
            disp(['***** Running analysis on: ' splitdir ' *****']);
            fNames = read_folder_contents(splitdir,'tif');

            [ scalinginfo, ~, lut ]=determine_scaling(splitdir, fNames, fullfile(lutfolder,lutfname) ,'degrees');

            Foveated_Montage_DFT_Analysis(splitdir, fNames, scalinginfo, 'degrees', lut, true, confocal_coords);
            something=true;
        end
                
        
        if ~something
            disp(['***** ' folderList{f} ' doesn''t contain any analyzable materials. *****']);
        end
    end
end

%% Process each of the above.
restartf=1;
endf=length(folderList);
%%
for f=restartf:endf
    restartf=f;
    confocal_coords=[];
    something=false;
    if isdir(f)
        confocaldir=fullfile(thisfolder, folderList{f},'confocal','Results_foveated');
        splitdir=fullfile(thisfolder, folderList{f},'split detection','Results_foveated');
        
        if exist(confocaldir, 'dir') && exist(splitdir, 'dir')
            fNameC = read_folder_contents(confocaldir,'mat');
            fNameS = read_folder_contents(splitdir,'mat');
            
            if ~isempty(fNameC) && ~isempty(fNameS)
                something=true;
                disp(['***** Merging data from: ' confocaldir ' and\n ' splitdir ' *****']);
                
                load(fullfile(confocaldir, fNameC{1}), 'density_map', 'blendederrim','fovea_coords');
                blendederrim_conf = blendederrim;
                density_map_conf = density_map;
                
                load(fullfile(splitdir, fNameS{1}), 'density_map', 'blendederrim');
                blendederrim_split = blendederrim;
                density_map_split = density_map;
                
                disp(['Data loaded, merging...']);
                clear density_map blendederrim
            
                errconfpolar = imcart2pseudopolar(blendederrim_conf, 1, .5, fovea_coords,'makima' , 0);
                errsplitpolar = imcart2pseudopolar(blendederrim_split, 1, .5, fovea_coords,'makima' , 0);
                errconfpolar(errconfpolar==0) = NaN;
                errsplitpolar(errsplitpolar==0) = NaN;
                
                figure(1); plot(mean(errconfpolar,'omitnan')); hold on; plot(mean(errsplitpolar,'omitnan')); hold off;
                drawnow;
                offset = 200;
                errconfpolar = errconfpolar(:, offset:end);
                errsplitpolar = errsplitpolar(:, offset:end);
                
                confannuli = zeros(size(blendederrim_conf));
                splitannuli = zeros(size(blendederrim_conf));

                avgdifferr = (mean(errconfpolar,'omitnan')-mean(errsplitpolar,'omitnan'))./ ...
                             ( (mean(errsplitpolar,'omitnan') + mean(errconfpolar,'omitnan'))/2 );

                         
                         
                confhigh  = find((avgdifferr>=-0.05) == 0, 1, 'first') + offset;
                splithigh = find( (avgdifferr<=0.05) == 1, 1, 'first') + offset; %Problem is here.

                blendrange = round((confhigh-splithigh)/2);
                confhigh = confhigh-blendrange;
                splithigh = splithigh+blendrange;

                confdisk = strel('disk',confhigh,0);
                confdisk = confdisk.Neighborhood;

                diskshiftx = floor(fovea_coords(1)-(size(confdisk,2)/2));
                diskshifty = floor(fovea_coords(2)-(size(confdisk,1)/2));

                confannuli(diskshifty:diskshifty+size(confdisk,1)-1,...
                           diskshiftx:diskshiftx+size(confdisk,2)-1) = confdisk;
                splitdisk = strel('disk',splithigh,0);
                splitdisk = splitdisk.Neighborhood;

                diskshiftx = floor(fovea_coords(1)-(size(splitdisk,2)/2));
                diskshifty = floor(fovea_coords(2)-(size(splitdisk,1)/2));

                splitannuli(diskshifty:diskshifty+size(splitdisk,1)-1,...
                           diskshiftx:diskshiftx+size(splitdisk,2)-1) = splitdisk;

                splitannuli = abs(1-splitannuli);

                confannuli = imgaussfilt(confannuli, blendrange);
                splitannuli = imgaussfilt(splitannuli, blendrange);

                % Weight each map against its filtered annuli
                density_map_conf = confannuli.*density_map_conf;
                density_map_split = splitannuli.*density_map_split;

                figure(1); imagesc(density_map_conf);
                figure(2); imagesc(density_map_split);

                density_map_comb = (density_map_conf+density_map_split)./(confannuli+splitannuli);
                
                denspolar = imcart2pseudopolar(density_map_comb, 1, .5, fovea_coords,'makima' , 0);
                denspolar(denspolar==0) = NaN;

                
                figure(3); imagesc(density_map_comb);
                drawnow;
                saveas(gcf, fullfile(thisfolder, folderList{f},'merged_density.png') );
                figure(4); plot(mean(denspolar, 'omitnan')); title('Merged density average');
                drawnow;
                saveas(gcf, fullfile(thisfolder, folderList{f},'merged_density_plot.png') );
            
            end
        end
        

        if ~something
            disp(['***** ' folderList{f} ' doesn''t contain any analyzable materials. *****']);
        end
    end
end

