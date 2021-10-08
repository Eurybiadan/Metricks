function [scaling_information, selectedunit, lut] = determine_scaling(basepath, fNames, lut_path, desired_unit)
% [scaling_information, selectedunit] = determine_scaling(basepath, fNames)
%   
% This function is responsible for determining the scale information for a
% future dataset.
% 
% It can take in a base path to look for a LUT file, and it can take in a
% list of filenames and determine their scaling information.
%
% Inputs:
%   @basepath: The path to start at when asking for a LUT file.
%
%   @fNames: A cell list of filenames that should be checked
%
% Outputs:
%   @scaling_information: If there are no filenames supplied, then this
%   will contain either a single double value corresponding to the user's
%   input scale, or the data from the selected LUT file.
%                         If there ARE filenames supplied, then this will
%   contain a list (equivalent in length to fNames) of scale values that
%   correspond to the same index in the fName list.
%
%   @selectedunit: The unit that the user decided they wanted.
%
%   @lut: The lut that was used to determine the scaling and selectedunit
%   information.
%
%   Created by Robert F Cooper, 2018-11-08

    if ~exist('basepath','var')
        basepath=pwd;
    end

    if ~exist('desired_unit','var')
        liststr = {'microns (mm density)','degrees','arcmin'};
        [selectedunit, oked] = listdlg('PromptString','Select output units:',...
                                      'SelectionMode','single',...
                                      'ListString',liststr);
        if oked == 0
            error('Cancelled by user.');
        end

        selectedunit = liststr{selectedunit};
    else
        selectedunit = desired_unit;
    end

    
    if ~exist('lut_path','var')
        [scaling_path, scalingdir] = uigetfile(fullfile(basepath,'*.csv'),'Select scaling LUT, OR cancel if you want to input the scale directly.');
        scaling_path = fullfile(scalingdir,scaling_path);
    else
        scaling_path = lut_path;
        
    end
    
    scaling_information = NaN;
    if scaling_path == 0        

        while isnan(scaling_information)                

            scaling_information = inputdlg('Input the scale in UNITS/PIXEL:','Input the scale in UNITS/PIXEL:');

            scaling_information = str2double(scaling_information);
            lut=[];
            if isempty(scaling_information)
                error('Cancelled by user.');
            end
        end
    else
        
        [~, lutData] = load_scale_file(scaling_path);

        lut= lutData;
        if ~exist('fNames','var')
            scaling_information = lutData;            
        else
            scaling_information = nan(length(fNames),1);
            
            for i=1:length(fNames)
                % Calculate the scale for this identifier.                                
                LUTindex=find( cellfun(@(s) ~isempty(strfind(fNames{i},s )), lutData{1} ) );

                % Use whichever scale is most similar to our filename.
                sim = 1000*ones(length(LUTindex),1);
                for l=1:length(LUTindex)
                    sim(l) = lev(fNames{i}, lutData{1}{LUTindex(l)});
                end
                [~,simind]=min(sim);
                LUTindex = LUTindex(simind);

                axiallength = lutData{2}(LUTindex);
                pixelsperdegree = lutData{3}(LUTindex);

                micronsperdegree = (291*axiallength)/24;

                switch selectedunit
                    case 'microns (mm density)'
                        scaling_information(i) = 1 / (pixelsperdegree / micronsperdegree);
                    case 'degrees'
                        scaling_information(i) = 1/pixelsperdegree;
                    case 'arcmin'
                        scaling_information(i) = 60/pixelsperdegree;
                end
            end
        end
    end
end

function [ scaling_row, scaling_col ] = load_scale_file( fileloc )
% Robert Cooper 06-18-2012
%   This function loads the needed scaling information for our metrics.

    fid = fopen(fileloc,'r');

    scaling_col = textscan(fid , '%s %f %f', 'delimiter', ',');

    scaling_col{1} = cellfun(@(s) strrep(s,'"',''),scaling_col{1},'UniformOutput',false);
    
    % Reform into cells that contain each row
    for i=1:size(scaling_col{1},1)


        ID = scaling_col{1}(i);
        axial{1} = scaling_col{2}(i);
        pix_per_deg{1} = scaling_col{3}(i);    

        scaling_row{i} = [ID axial pix_per_deg];

    end

    fclose(fid);
end

