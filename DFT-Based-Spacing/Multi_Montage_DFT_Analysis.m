% Robert F Cooper
% 1/13/21
%
% This script's sole job is to run multiple montages in series; to save time, 
% instead of operating on individual files though, it operates on a series
% of folders.

folderList = {}

thisfolder = pwd;


while thisfolder ~= 0

    thisfolder = uigetdir(thisfolder, 'Select the folders containing the montages you wish to assess.');
    
    if thisfolder ~= 0
        folderList = [folderList; thisfolder];
    end
end

[lutfname, lutfolder] = uigetfile(fullfile(pwd,'*.csv'),'Select scaling LUT, OR cancel if you want to input the scale directly.');

%%
for f=1:length(folderList)
    
    fNames = read_folder_contents(folderList{f},'tif');
    
    [ scalinginfo, ~, lut ]=determine_scaling(folderList{f}, fNames, fullfile(lutfolder,lutfname) ,'degrees');
    
    
    Montage_DFT_Analysis(folderList{f}, fNames, scalinginfo, 'degrees', lut, true);
    
end