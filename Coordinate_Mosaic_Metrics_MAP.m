% Robert Cooper
% 2017-09-29
%
% This script calculates the coordinate metrics from a selected folder.


clear;
close all force;

WINDOW_SIZE = 35;

%% Crop the coordinates/image to this size in [scale], and calculate the area from it.
% If left empty, it uses the size of the image.

if length(WINDOW_SIZE) > 1 && ~isempty(WINDOW_SIZE)
   error('Window size cannot be empty, and must be a SINGLE value!');
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

            pixelwindowsize = WINDOW_SIZE/scaleval;
            
            if exist(fullfile(basepath, [fnamelist{i}(1:end-length('_coords.csv')) '.tif']), 'file')

                im = imread( fullfile(basepath, [fnamelist{i}(1:end-length('_coords.csv')) '.tif']));

                width = size(im,2);
                height = size(im,1);
                maxrowval = height;
                maxcolval = width;
            else

                width  = max(coords(:,1)) - min(coords(:,1));
                height = max(coords(:,2)) - min(coords(:,2));
                maxrowval = max(coords(:,2));
                maxcolval = max(coords(:,1));
            end
            
            statistics = cell(size(coords,1),1);
            
            tic;
            parfor c=1:size(coords,1)
                
                rowborders = round([coords(c,2)-(pixelwindowsize/2) coords(c,2)+(pixelwindowsize/2)]);
                colborders = round([coords(c,1)-(pixelwindowsize/2) coords(c,1)+(pixelwindowsize/2)]);

                rowborders(rowborders<1) =1;
                colborders(colborders<1) =1;
                rowborders(rowborders>maxrowval) =maxrowval;
                colborders(colborders>maxcolval) =maxcolval;
                
                clipped_coords =coordclip(coords,colborders,...
                                                 rowborders,'i');
                % [xmin xmax ymin ymax] 
                clip_start_end = [colborders rowborders];

                statistics{c} = determine_mosaic_stats( clipped_coords, scaleval, selectedunit, clip_start_end ,[colborders(2)-colborders(1) rowborders(2)-rowborders(1)], 4 );

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Determine FFT Power Spectra %%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if (exist('fit_fourier_spacing.m','file') == 2) && exist(fullfile(basepath, [fnamelist{i}(1:end-length('_coords.csv')) '.tif']), 'file')==2
                    [pixel_spac, interped_map] = fit_fourier_spacing(im);
                    statistics{c}.DFT_Spacing = pixel_spac*scaleval;                
                end


                warning off;
                [ success ] = mkdir(basepath,'Results');
                warning on;
            
            end
            toc;
            
            %% Map output
            metriclist = fieldnames(statistics{1});
            [selectedmetric, oked] = listdlg('PromptString','Select map metric:',...
                                          'SelectionMode','single',...
                                          'ListString',metriclist);
            
            if oked == 0
                error('Cancelled by user.');
            end
                                      
            interped_map=zeros([height width]);
            sum_map=zeros([height width]);

            
            for c=1:size(coords,1)

                    thisval = statistics{c}.(metriclist{selectedmetric}); 

                    rowrange = round(coords(c,2)-(pixelwindowsize/2):coords(c,2)+(pixelwindowsize/2));
                    colrange = round(coords(c,1)-(pixelwindowsize/2):coords(c,1)+(pixelwindowsize/2));

                    rowrange(rowrange<1) =[];
                    colrange(colrange<1) =[];
                    rowrange(rowrange>maxrowval) =[];
                    colrange(colrange>maxcolval) =[];
                    
                    interped_map(rowrange,colrange) = interped_map(rowrange,colrange) + thisval;
                    sum_map(rowrange, colrange) = sum_map(rowrange, colrange) + 1;

            end
                %

            interped_map = interped_map./sum_map;

            dispfig=figure(1); imagesc(interped_map); axis image; colorbar;
            title(['Minimum value: ' num2str(min(interped_map(:))) ' Maximum value: ' num2str(max(interped_map(:)))])
            
            result_fname = [getparent(basepath,'short') '_density_bound_coordmap_' date '_' num2str(WINDOW_SIZE) '_' metriclist{selectedmetric}];
            
            saveas(gcf,fullfile(basepath,'Results', [result_fname '_fig.png']));
            saveas(gcf,fullfile(basepath,'Results', [result_fname '_fig.svg']));
            
            scaled_map = interped_map-min(interped_map(:));
            scaled_map = uint8(255*scaled_map./max(scaled_map(:)));
            imwrite(scaled_map, parula(256), fullfile(basepath,'Results', [result_fname '_raw.tif']))

            %%
        end
    catch ex
        warning(['Unable to analyze ' fnamelist{i} ':']);
        warning([ex.message ', In file: ' ex.stack(1).file '  Line: ' num2str(ex.stack(1).line)]);
    end
end
close(proghand);
