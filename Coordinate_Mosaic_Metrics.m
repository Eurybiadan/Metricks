% Copyright (C) 2019 Robert F Cooper
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% Metricks - A MATLAB package for analyzing the cone photoreceptor mosaic.
% 
% Coordinate_Mosiac_Metrics calculates the metrics for every
% image/coordinate pair in a given folder.
%
% When run, the script will prompt the user to select a folder with image/coordinate pairs.
% 
% **At present, images must be 8-bit grayscale tifs, coordinates must be formatted 
%   as a 2 column matrix (x,y), and must be named using the following convention,
%   where [imagename] can be any valid filename:**
% * Image File: [imagename].tif
% * Coordinate File: [imagename]\_coords.csv
% 
% It will then prompt the user to select what the output unit should be. At present,
% the options are:
% * Microns (using millimeters^2 for density)
% * Degrees
% * Arcminutes
% 
% Once the output unit is select, it will give the user the option to pick a 
% lookup table. The lookup table allows the software to analyze a folder of 
% images from different subjects/timepoints/conditions. The lookup table itself
% **must** be a 3 column 'csv' file, where the **first column** is a common 
% identifier for image/coordinate pairs, the **second column** is the axial 
% length (or '24' if the axial length is unknown) of the image/coordinate pairs,
% and the **third column** is the pixels per degree of the image/coordinate pairs.
% Each row must contain a different identifier/axial length/pixels per degree tuple.
% 
% An example common identifier could be a subject number, e.g, when working with the files:
% - 1235_dateoftheyear_OD_0004.tif
% - 1235_dateoftheyear_OD_0005.tif
% 
% Common identifiers could be "1235", "1235_dateoftheyear", "1235_dateoftheyear_OD".
% If all three were placed in a LUT, then the one that matches the most (as determined
% via levenshtein distance) will be used. In this example, we would use "1235_dateoftheyear_OD".
% 
% If we had another date, say: 1235_differentdateoftheyear_OD_0005.tif, then 
% _only_ the identifier "1235" would match between all images. However, say the
% two dates have different scales, then you would want to create two rows in the
% look up table for each date, with identifiers like: "1235_dateoftheyear" and
% "1235_differentdateoftheyear".
% 
% **If you do not wish to use a lookup table, then press "cancel", and the software
% will allow you put in your own scale in UNITS/pixel.**
% 
% **This software has the ability to pre-crop the input data (if, for example,
% you have 80 pixels of coordinates and you only want to analyze the middle 50).**
% 
% To specify a cropping window, input the size (in the units you are going to 
% use) in to the brackets on line ~129 (of the variable 
% windowsize) of Coordinate_Mosaic_Metrics.m.
% 
% **Cropping is governed by the following rules:**
% 
% 1) If the tif is present and windowsize is not specified, the analysis will 
% be done on everything within the dimensions of the image.
% 2) If the tif is present and windowsize is specified, the assumed center of 
% the image is calculated according to the borders of the tif. **In either case,
% it doesn’t “care” how many (or even if there are any) cells in the image.**
% 3) If the tif is not present and windowsize is not specified, the analysis will
% be done on everything within the min and max coordinates in both x and y directions.
% So if you have an image in which there is an absence of cells on one side, 
% for example, you might end up with a clipped area that is not a square.
% 4) If the tif is not present and windowsize is specified, the assumed center 
% of the image is calculated according to the min and max coordinates in both 
% x and y directions. So if you have an image in which there is an absence of 
% cells on one side, the center will shift towards the other side of the image.
% 
% 
% The software will then run, and calculate every metric currently validated.
% 
% At present, it calculates the following metrics from each image and coordinate pair:
% 
% - Number of Unbound Cells
% - Number of Bound Cells
% - Total Area
% - Total Bounded Area
% - Mean Voronoi Area
% - Percent Six-Sided Voronoi
% - Density (uncorrected/corrected)
% - Nearest Neighbor Distance (uncorrected/corrected)
% - Inter-Cell Distance (uncorrected/corrected)
% - Furthest Neighbor Distance (uncorrected/corrected)
% - Density Recovery Profile Distance
% - Voronoi Area Regularity Index
% - Voronoi Number of Sides Regularity Index
% - Nearest Neighbor Regularity Index
% - Inter-Cell Regularity Index
% 
% The results will then be placed in to a datestamped file within a "Results" 
% folder as a subfolder of the one selected for analysis.
% 
%
% Don't thank me; cite me:
% 
% Every metric that is run via the main "Coordinate_Mosaic_Metrics.m" script
% has been validated and used in the following manuscript: 
% 
% Cooper RF, Wilk MA, Tarima S, Dubra A, Carroll J. 
% “Evaluating descriptive metrics of the human cone mosaic.”
% Invest Ophthalmol Vis Sci. 2016 57(7):2993.
% 
% You can also find formal definitions of each metric calculated here in that paper.
% 
% **This package is free for use under GPL v3, but I ask that you please cite 
% the above paper if you use this package.**
%
% 


clear;
close all force;

windowsize = [];
%% Crop the coordinates/image to this size in [scale], and calculate the area from it.
% If left empty, it uses the size of the image.

if length(windowsize) > 1
   error('Window size can only be empty ([]), or a single value!');
end

basePath = which('Coordinate_Mosaic_Metrics.m');

[basePath ] = fileparts(basePath);
path(path,fullfile(basePath,'lib')); % Add our support library to the path.

[basepath] = uigetdir(pwd);

[fnamelist, isdir ] = read_folder_contents(basepath,'csv');
[fnamelisttxt, isdirtxt ] = read_folder_contents(basepath,'txt');

fnamelist = [fnamelist; fnamelisttxt];
isdir = [isdir;isdirtxt];

liststr = {'microns (mm density)','degrees','arcmin'};
[selectedunit, oked] = listdlg('PromptString','Select output units:',...
                              'SelectionMode','single',...
                              'ListString',liststr);
if oked == 0
    error('Cancelled by user.');
end

selectedunit = liststr{selectedunit};                          

[scalingfname, scalingpath] = uigetfile(fullfile(basepath,'*.csv'),'Select scaling LUT, OR cancel if you want to input the scale directly.');

scaleinput = NaN;
if scalingfname == 0        
    
    while isnan(scaleinput)                
        
        scaleinput = inputdlg('Input the scale in UNITS/PIXEL:','Input the scale in UNITS/PIXEL:');
        
        scaleinput = str2double(scaleinput);
        
        if isempty(scaleinput)
            error('Cancelled by user.');
        end
    end
else
    [~, lutData] = load_scaling_file(fullfile(scalingpath,scalingfname));
end


first = true;

proghand = waitbar(0,'Processing...');

for i=1:size(fnamelist,1)

    try
        if ~isdir{i}

            
            if length(fnamelist{i})>42
                waitbar(i/size(fnamelist,1), proghand, strrep(fnamelist{i}(1:42),'_','\_') );
            else
                waitbar(i/size(fnamelist,1), proghand, strrep(fnamelist{i},'_','\_') );
            end

            if isnan(scaleinput)
                % Calculate the scale for this identifier.                                
                LUTindex=find( cellfun(@(s) ~isempty(strfind(fnamelist{i},s )), lutData{1} ) );

                % Use whichever scale is most similar to our filename.
                sim = 1000*ones(length(LUTindex),1);
                for l=1:length(LUTindex)
                    sim(l) = lev(fnamelist{i}, lutData{1}{LUTindex(l)});
                end
                [~,simind]=min(sim);
                LUTindex = LUTindex(simind);
                
                axiallength = lutData{2}(LUTindex);
                pixelsperdegree = lutData{3}(LUTindex);

                micronsperdegree = (291*axiallength)/24;
                
                switch selectedunit
                    case 'microns (mm density)'
                        scaleval = 1 / (pixelsperdegree / micronsperdegree);
                    case 'degrees'
                        scaleval = 1/pixelsperdegree;
                    case 'arcmin'
                        scaleval = 60/pixelsperdegree;
                end
            else
                scaleval = scaleinput;
            end


            %Read in coordinates - assumes x,y
            coords=dlmread(fullfile(basepath,fnamelist{i}));
            
            % It should ONLY be a coordinate list, that means x,y, and
            % nothing else.
            if size(coords,2) ~= 2
                warning('Coordinate list contains more than 2 columns! Skipping...');
                continue;
            end

            % If the corresponding image exists in the folder, use the image bounds to calculate our sizes
            if exist(fullfile(basepath, [fnamelist{i}(1:end-length('_coords.csv')) '.tif']), 'file')

                im = imread( fullfile(basepath, [fnamelist{i}(1:end-length('_coords.csv')) '.tif']));

                width = size(im,2);
                height = size(im,1);

                if ~isempty(windowsize)
                    pixelwindowsize = windowsize/scaleval;

                    diffwidth  = (width-pixelwindowsize)/2;
                    diffheight = (height-pixelwindowsize)/2;
                    
                    if diffwidth<0
                        diffwidth=0;
                    end                    
                    if diffheight<0
                        diffheight=0;
                    end
                else

                    pixelwindowsize = [height width];
                    diffwidth=0;
                    diffheight=0;
                end

                clipped_coords =coordclip(coords,[diffwidth  width-diffwidth],...
                                                 [diffheight height-diffheight],'i');

                clip_start_end = [diffwidth+1  width-diffwidth diffheight+1 height-diffheight];
            else

                width  = max(coords(:,1)) - min(coords(:,1));
                height = max(coords(:,2)) - min(coords(:,2));

                if ~isempty(windowsize)
                    pixelwindowsize = windowsize/scaleval;

                    diffwidth  = (width-pixelwindowsize)/2;
                    diffheight = (height-pixelwindowsize)/2;
                else
                    pixelwindowsize = [height width];
                    diffwidth=0;
                    diffheight=0;
                end

                clipped_coords =coordclip(coords,[min(coords(:,1))+diffwidth  max(coords(:,1))-diffwidth],...
                                                 [min(coords(:,2))+diffheight max(coords(:,2))-diffheight],'i');

                clip_start_end = [min(coords(:,1))+diffwidth  max(coords(:,1))-diffwidth min(coords(:,2))+diffheight max(coords(:,2))-diffheight];
            end


            statistics = determine_mosaic_stats( clipped_coords, scaleval, selectedunit, clip_start_end ,[pixelwindowsize pixelwindowsize], 4 );

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Determine FFT Power Spectra %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (exist('fit_fourier_spacing') == 2) && exist(fullfile(basepath, [fnamelist{i}(1:end-length('_coords.csv')) '.tif']), 'file')==2
                
                clipped_im = im(round(clip_start_end(3):clip_start_end(4)), round(clip_start_end(1):clip_start_end(2)) );
                
                [pixel_spac, ~, quality] = fit_fourier_spacing(clipped_im, min(size(clipped_im)), false,'row');
                statistics.DFT_Row_Spacing = pixel_spac*scaleval;
                statistics.DFT_Row_Quality = quality;
                
                [pixel_spac, ~, quality] = fit_fourier_spacing(clipped_im, min(size(clipped_im)), false,'cell');
                statistics.DFT_Cell_Spacing = pixel_spac*scaleval;
                statistics.DFT_Cell_Quality = quality;
                
                
            end


            warning off;
            [ success ] = mkdir(basepath,'Results');
            warning on;
            
            if isempty(windowsize)
                result_fname = [getparent(basepath,'short') '_coordstats_' date '.csv'];
            else
                result_fname = [getparent(basepath,'short') '_coordstats_' date '_' num2str(windowsize) selectedunit '.csv'];
            end
            if success

                if first
                    fid= fopen(fullfile(basepath,'Results', result_fname),'w');

                    % If it is the first time writing the file, then write the
                    % header
                    fprintf(fid,'Filename');

                    % Grab the names of the fields we're working with
                    datafields = fieldnames(statistics);

                    numfields = size(datafields,1);                

                    k=1;

                    while k <= numfields

                        val = statistics.(datafields{k});

                        % If it is a multi-dimensional field, remove it
                        % from our csv, and write it separately.
                        if size(val,1) ~= 1 || size(val,2) ~= 1   
                            disp([datafields{k} ' removed!']);
                            datafields = datafields([1:k-1 k+1:end]);                        
                            numfields = numfields-1;                        
                        else
    %                         disp([fields{k} ' added!']);
                            fprintf(fid,',%s',datafields{k});
                            k = k+1;
                        end 


                    end  
                    fprintf(fid,'\n');

                    first = false;

                else % If it isn't the first entry, then append.
                    fid= fopen(fullfile(basepath,'Results',result_fname ),'a');
                end

                % Write the file we've worked on as the first column
                fprintf(fid,'%s', fnamelist{i});

                for k=1:size(datafields,1)
    %                 fields{k}
                    if size(val,1) == 1 || size(val,2) == 1
                        val = statistics.(datafields{k});

                        fprintf(fid,',%1.2f',val);
                    end
                end

                fprintf(fid,'\n');
                fclose(fid);
            else
                error('Failed to make results folder! Exiting...');
            end

        end
    catch ex
        warning(['Unable to analyze ' fnamelist{i} ':']);
        warning([ex.message ', In file: ' ex.stack(1).file '  Line: ' num2str(ex.stack(1).line)]);
    end
end
close(proghand);
