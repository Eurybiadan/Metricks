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
restartf=21;
endf=21;
% endf=length(folderList);
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