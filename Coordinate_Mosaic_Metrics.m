% Robert Cooper
% 2017-09-29
%
% This script calculates the coordinate metrics from a selected folder.


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

            waitbar(i/size(fnamelist,1), proghand, strrep(fnamelist{i}(1:42),'_','\_') );

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
                else

                    pixelwindowsize = [height width];
                    diffwidth=0;
                    diffheight=0;
                end

                clipped_coords =coordclip(coords,[diffwidth  width-diffwidth],...
                                                 [diffheight height-diffheight],'i');

                clip_start_end = [diffheight height-diffheight diffwidth  width-diffwidth];
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

                clip_start_end = [min(coords(:,2))+diffheight max(coords(:,2))-diffheight min(coords(:,1))+diffwidth  max(coords(:,1))-diffwidth];
            end


            statistics = determine_mosaic_stats( clipped_coords, scaleval, selectedunit, clip_start_end ,[pixelwindowsize pixelwindowsize], 4 );

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Determine FFT Power Spectra %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (exist('fit_fourier_spacing') == 2) && exist(fullfile(basepath, [fnamelist{i}(1:end-length('_coords.csv')) '.tif']), 'file')==2
                [pixel_spac, interped_spac_map] = fit_fourier_spacing(im);
                statistics.DFT_Spacing = pixel_spac*scaleval;                
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
