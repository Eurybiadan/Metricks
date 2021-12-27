%% Robert F Cooper
% 1/13/21
%
% This script's sole job is to run multiple montages in series; to save time, 
% instead of operating on individual files though, it operates on a series
% of folders.
clear;
close all;
folderList = {};

thisfolder = pwd;



thisfolder = uigetdir(thisfolder, 'Select the folders containing the montages you wish to assess.');
    
[folderList, isdir]=read_folder_contents(thisfolder);

[lutfname, lutfolder] = uigetfile(fullfile(pwd,'*.csv'),'Select scaling LUT, OR cancel if you want to input the scale directly.');

%%
restartf = 45;
endf = 45; %length(folderList);
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

            [confocal_coords, mask]=Foveated_Montage_DFT_Analysis(confocaldir, fNames, scalinginfo, 'degrees', lut, true);
            something=true;
        end
        close all;
        
        if exist(splitdir, 'dir')
            disp(['***** Running analysis on: ' splitdir ' *****']);
            fNames = read_folder_contents(splitdir,'tif');

            [ scalinginfo, ~, lut ]=determine_scaling(splitdir, fNames, fullfile(lutfolder,lutfname) ,'degrees');

            Foveated_Montage_DFT_Analysis(splitdir, fNames, scalinginfo, 'degrees', lut, true, confocal_coords, mask);
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
                
                load(fullfile(confocaldir, fNameC{1}), 'density_map', 'blendederrim','fovea_coords', 'scaling');
                blendederrim_conf = blendederrim;
                density_map_conf = density_map;
                confocalmask = ~isnan(density_map_conf);
                
                load(fullfile(splitdir, fNameS{1}), 'density_map', 'blendederrim');
                blendederrim_split = blendederrim.*confocalmask;
                density_map_split = density_map.*confocalmask;
                
                disp(['Data loaded, merging...']);
                clear density_map blendederrim
            
                errconfpolar = imcart2pseudopolar(blendederrim_conf, 1, .5, fovea_coords,'makima' , 0);
                errsplitpolar = imcart2pseudopolar(blendederrim_split, 1, .5, fovea_coords,'makima' , 0);
                errconfpolar(errconfpolar<=0) = NaN;
                errsplitpolar(errsplitpolar<=0) = NaN;
                
                
                drawnow;
                
                
                confannuli = zeros(size(blendederrim_conf));
                splitannuli = zeros(size(blendederrim_conf));

                avgdifferr = (mean(errconfpolar,'omitnan')-mean(errsplitpolar,'omitnan'))./ ...
                             ( (mean(errsplitpolar,'omitnan') + mean(errconfpolar,'omitnan'))/2 );

                figure(5); plot(mean(errconfpolar,'omitnan')); hold on; plot(mean(errsplitpolar,'omitnan')); plot(avgdifferr); hold off;         
                axis([0 length(errconfpolar) -2 2]); legend('Confocal Confidence', 'Split Confidence', '% diff');
                drawnow;
                saveas(gcf, fullfile(thisfolder, folderList{f},'conf_split_polarerror.png') );
                
                offset = find((avgdifferr>=0) == 1, 1, 'first');    
                
                MAXOFFSET = 2000;
                
                % We want to make sure that we found a valid blending spot.
                if ~isempty(offset) && offset < MAXOFFSET
                    
                    confhigh  = find((avgdifferr(:, offset:end)>=-0.2) == 0, 1, 'first') + offset;
                    splithigh = find( (avgdifferr(:, offset:end)<=0.2) == 1, 1, 'first') + offset;

                    if isempty(confhigh)
                        confhigh = find((avgdifferr(:, offset:end)>=0) == 0, 1, 'first') + offset;
                    end
                    
                    if isempty(splithigh)
                        splithigh = find((avgdifferr(:, offset:end)<=0) == 0, 1, 'first') + offset;
                    end
                    
                    blendrange = round((confhigh-splithigh)/2);
                                                    
                    mergeloc = confhigh-blendrange;

                else
                    disp(['Warning: unable to find ideal merging location. Guessing from closest values...']);
                    [val, offset] = max(avgdifferr(1:MAXOFFSET));
                    confhigh  = offset; %this needs work.
                    mergeloc = confhigh;
                    blendrange = 512;
                end
                

                confdisk = strel('disk',mergeloc,0);
                confdisk = confdisk.Neighborhood;

                diskshiftx = floor(fovea_coords(1)-(size(confdisk,2)/2));
                disksizex = diskshiftx+size(confdisk,2)-1;
                diskshifty = floor(fovea_coords(2)-(size(confdisk,1)/2));
                disksizey = diskshifty+size(confdisk,1)-1;
                
                diskstartx = 1;
                diskstarty = 1;
                diffx=0;    
                diffy=0;
                
                if diskshiftx<1
                    diskstartx = (-diskshiftx)+1;
                    diskshiftx=1;                    
                end
                if diskshifty<1
                    diskstarty = (-diskshifty)+1;
                    diskshifty=1;                    
                end
                
                if disksizex > size(confannuli,2)
                    diffx = size(confannuli,2)-disksizex;
                    disksizex = disksizex+diffx;
                end
                if disksizey > size(confannuli,1)
                    diffy = size(confannuli,1)-disksizey;
                    disksizey = disksizey+diffy;
                end
                
                confannuli(diskshifty:disksizey,...
                           diskshiftx:disksizex) = confdisk(diskstarty:end+diffy,...
                                                            diskstartx:end+diffx);
                       
                splitannuli = abs(1-confannuli);

                figure(1); subplot(1,2,1); imagesc(confannuli); axis image;
                
                confannuli = imgaussfilt(confannuli, blendrange/2);
                splitannuli = imgaussfilt(splitannuli, blendrange/2);

                subplot(1,2,2); imagesc(confannuli); axis image;
                
                % Weight each map against its filtered annuli
                density_map_conf = confannuli.*density_map_conf;
                density_map_split = splitannuli.*density_map_split;

                
                figure(2); subplot(1,2,1); imagesc(density_map_split); axis image;
                subplot(1,2,2); imagesc(density_map_conf);axis image;
                drawnow;
                saveas(gcf, fullfile(thisfolder, folderList{f},'conf_split_density_areas.png') );
                
                
                density_map_comb = (density_map_conf+density_map_split)./(confannuli+splitannuli);
                
                clear density_map_conf density_map_split
                
%                 denspolar = imcart2pseudopolar(density_map_comb, 1, .5, fovea_coords,'makima' , 0);
%                 denspolar(denspolar==0) = NaN;

                blendederrim_conf = confannuli.*blendederrim_conf;
                blendederrim_split = splitannuli.*blendederrim_split;
                
                blendederr_comb = (blendederrim_conf+blendederrim_split)./(confannuli+splitannuli);
                
                % Find our new lowest low, ignoring NaNs.
%                 highconf = quantile(blendederr_comb(~isnan(blendederr_comb)), 0.01);
%                 
%                 density_map_comb(blendederr_comb<=highconf) = NaN;     
%                 blendederr_comb(blendederr_comb<=highconf) = NaN;
                
                
                [X, Y] = meshgrid(1:size(density_map_comb,2), 1:size(density_map_comb,1));
                
                % Find all the points within 1.75 degrees of the fovea
                inrangemask = sqrt((X-fovea_coords(1)).^2 + (Y-fovea_coords(2)).^2) <= 1.75/scaling;
                % Remove all values outside the fovea that are more than 66% of the maximal foveal values.
                
                nohigher = max(max(density_map_comb.*inrangemask))*.66
                
                density_map_comb(( (density_map_comb > nohigher) & ~inrangemask )) = NaN;
                           
                
                clear blendederrim_split blendederrim_conf 
                
                
                figure(3); imagesc(density_map_comb); axis image;
                drawnow;
%                 pause;
                saveas(gcf, fullfile(thisfolder, folderList{f},'merged_density.png') );
%                 figure(4); plot(mean(denspolar, 'omitnan')); title('Merged density average');
                imagesc(blendederr_comb); axis image;
                drawnow;
                saveas(gcf, fullfile(thisfolder, folderList{f},'merged_error.png') );
%                 close all;
                
                
                safesave(fullfile(confocaldir, fNameC{1}), fullfile(thisfolder, folderList{f}, strrep(fNameC{1}, 'confocal', 'merged')), density_map_comb, blendederr_comb, confannuli, splitannuli);
                safesave_sm(fullfile(confocaldir, fNameC{1}), fullfile(thisfolder, 'Aggregation_Analysis', strrep(fNameC{1}, 'confocal', 'merged')), density_map_comb, blendederr_comb);
                clear density_map_comb confannuli splitannuli
            end
        end
        

        if ~something
            disp(['***** ' folderList{f} ' doesn''t contain any analyzable materials. *****']);
        end
    end
end

function []= safesave_sm(baseload, newsave, denscomb, errcomb)
    disp(['Using path: ' baseload ' as base for ' newsave]);
    load(baseload);
    density_map_comb = denscomb;
    blendederr_comb = errcomb;
    save(newsave, 'density_map_comb', 'blendederr_comb', 'fovea_coords', 'imsize', 'scaling',  '-v7.3')
end

function []= safesave(baseload, newsave, denscomb, errcomb, confannuli, splitannuli)
    disp(['Using path: ' baseload ' as base for ' newsave]);
    load(baseload);
    density_map_comb = denscomb;
    blendederr_comb = errcomb;
    confocal_annulus = confannuli;
    split_annulus = splitannuli;
    save(newsave, '-v7.3')
end

