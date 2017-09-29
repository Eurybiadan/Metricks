% close all;


pathname = pwd;
[fNames] = read_folder_contents(pathname,'csv');

% Load the scaling file
[tmp, lutData] = load_scaling_file('R:\Rob Cooper\Cell_Orientation_2015\LUT.csv');
% [tmp, lutData] = load_scaling_file('M:\Documents\Grad school\Cell_Orientation_2015\LUT.csv');
% [tmp, lutData] = load_scaling_file('/remote_project_folders/Rob Cooper/Cell_Orientation_2015/LUT.csv');

numnames = length(fNames);
avgpatchsize=0;
for n=1:numnames
    n
    fname = fNames{n}
    imname = [fname(1:end-11) '.tif'];
    % Read in the image, if it exists.
    impath = fullfile(pathname,imname);
    if exist(impath,'file')
       im = imread(impath);
    end
    
    % Calculate scale for this file
    [idpiece1 remain]=strtok(fname,'_'); %Take Referrer
    [idpiece2 remain]=strtok(remain,'_'); %Take ID #
    subID=[idpiece1 '_' idpiece2]; 
    clear remain idpiece1 idpiece2;
    
    LUTindex=find(strcmp(lutData{1},subID));

    axiallength = lutData{2}(LUTindex);
    pixelsperdegree = lutData{3}(LUTindex);

    micronsperdegree = (291*axiallength)/24;
    um_per_pix = 1 / (pixelsperdegree / micronsperdegree);
%     um_per_pix = .45;

    coords = dlmread(fullfile(pathname,fname));    
    
    % Calculate the cell distances
    mean_nn_dist=calc_icd(coords,um_per_pix);

    if ~exist('im','var')
        [im, shiftcoords] = generate_sim_view( coords ); 
        im = uint8(im);

        coords = shiftcoords;

    end
%     if n==11
%        n 
%     end
%     [pumpatches{n}, pum_opatches{n}] = find_connected_patches( im, coords, PUM{n});
%     [rfcpatches{n}, rfc_opatches{n}] = find_connected_patches( im, coords, RFC{n});
%     [peppatches{n}, pep_opatches{n}] = find_connected_patches( im, coords, PEPPE{n});

    % Find the patch borders
    imsize = size(im);
    patchsize = imsize / 8;
    patchwidths = 0:patchsize:imsize;
    
    kk=1;
    for ii=1:(length(patchwidths)-1)
        beginpatchy = patchwidths(ii);
        endpatchy   = ceil(patchwidths(ii+1));
        for jj=1:(length(patchwidths)-1)
            beginpatchx = patchwidths(jj);
            endpatchx   = ceil(patchwidths(jj+1));
        
    
            [patchcoords{n}{kk} patches{n}{kk}] = coordclipv2(coords, [beginpatchx endpatchx], [beginpatchy endpatchy]);
            kk=kk+1;
        end
    end
    
    for ii=1:length(patches{n})
        avgpatchsize = avgpatchsize + length(patches{n}{ii});
        pum_opatches{n}{ii} = PUM{n}(patches{n}{ii});
        pum_opatches{n}{ii} = pum_opatches{n}{ii}(pum_opatches{n}{ii} ~= -100);
        
        rfc_opatches{n}{ii} = RFC{n}(patches{n}{ii});
        rfc_opatches{n}{ii} = rfc_opatches{n}{ii}(rfc_opatches{n}{ii} ~= -100);
        
        pep_opatches{n}{ii} = PEPPE{n}(patches{n}{ii});
        pep_opatches{n}{ii} = pep_opatches{n}{ii}(pep_opatches{n}{ii} ~= -100);
    end

%     for i=1:length(pumpatches{n})
%         rfc_opatches{n}{i} = RFC{n}(pumpatches{n}{i});    
%     end
%     
%     for i=1:length(pumpatches{n})
%         pep_opatches{n}{i} = PEPPE{n}(pumpatches{n}{i});    
%     end
    
    % Analyzing %6 sides
    triangulation = delaunayTriangulation(coords);
    numsixsided =0;
    for i=1:length(triangulation.Points)

        neigh = cell2mat( vertexAttachments(triangulation,i) );

        if size(neigh,2) == 6
            numsixsided = numsixsided +1;
        end

    end
    
    percentsixsided(n) = 100*numsixsided./length(coords);
    mean_icd(n) = mean_nn_dist;
    clear im;
   
end


toanalyze = percentsixsided > 65;

PatchPooledVarianceAnalysis

pumnumsmall
rfcnumsmall
pepnumsmall

%% Setting up for SPSS insertion
rfcall = [];
pumall = [];
pepall = [];
% Analyze their differences.
for i=1:length(fNames)

    
    if toanalyze(i)
        i
    % For comparing submosaics and real mosaics
    % Find all cases where they match
    matches = logical((RFC{i}~=-100) & (PUM{i}~=-100) & (PEPPE{i}~=-100));
    rfc = (RFC{i}(matches));
    pum = (PUM{i}(matches));
    pep = (PEPPE{i}(matches));
    
    
 
%     diffrfcpum = (rfc-pum);
%     diffrfcpep = (rfc-pep);
%     diffpeppum = (pep-pum);
%     
%     diffs = [diffrfcpum diffrfcpep diffpeppum];
%     phasewrapped = abs(diffs) > 30;
%     
%     % Determine our biggest difference- this will be our anchor
%     [vals inds]=max( abs( diffs ),[],2);
%     
%     % Shift the difference to the anchor
%     for j=1:length(inds)
%         
%         if phasewrapped(j, inds(j))
%         switch inds(j)
%             case 1 % between rfc and pum
%                 if sign(diffrfcpum(j)) == 1 % rfc is higher, and the anchor
%                     pum(j) = pum(j) +60;
%                     
%                     % Flip pep over the
%                     % anchor as well to maintain relative distances
%                     if abs(diffrfcpep(j)) > 30
%                         pep(j) = pep(j)+60;
%                     else                                                
%                         pep(j) = diffrfcpep(j).*2+pep(j);
%                     end
%                     
%                 else % pum higher, and is the anchor- bring rfc to it
%                     rfc(j) = rfc(j) +60;
%                     
%                     % Flip pep over the
%                     % anchor as well to maintain relative distances
%                     if abs(diffpeppum(j)) > 30
%                         pep(j) = pep(j)+60;
%                     else
%                         pep(j) = -diffpeppum(j).*2+pep(j);
%                     end
%                     
%                 end
% 
%             case 2 % between rfc and pep
%                 if sign(diffrfcpep(j)) == 1 % rfc is the anchor
%                     pep(j) = pep(j) +60;
%                                      
%                     % Flip pep over the
%                     % anchor as well to maintain relative distances
%                     if abs(diffrfcpum(j)) > 30
%                         pum(j) = pum(j)+60;
%                     else                      
%                         pum(j) = diffrfcpum(j).*2+pum(j);                        
%                     end                   
%                     
%                 else % pep is the anchor
%                     rfc(j) = rfc(j) +60;
%                     
%                     % Flip pep over the
%                     % anchor as well to maintain relative distances
%                     if abs(diffpeppum(j)) > 30
%                         pum(j) = pum(j)+60;
%                     else
%                         pum(j) = diffpeppum(j).*2+pum(j);
%                     end 
%                     
%                 end
%             case 3 % between pep and pum
%                 if sign(diffpeppum(j)) == 1 % pep is the anchor
%                     pum(j) = pum(j) +60;
%                                        
%                     % Flip pep over the
%                     % anchor as well to maintain relative distances
%                     if abs(diffrfcpep(j)) > 30
%                         rfc(j) = rfc(j)+60;
%                     else
%                         rfc(j) = -diffrfcpep(j).*2+rfc(j);
%                     end                   
%                     
%                 else % pum is the anchor
%                     pep(j) = pep(j) +60;
%                     
%                     % Flip pep over the
%                     % anchor as well to maintain relative distances
%                     if abs(diffrfcpum(j)) > 30
%                         rfc(j) = rfc(j)+60;
%                     else
%                         rfc(j) = -diffrfcpum(j).*2+rfc(j);
%                     end 
%                     
%                 end
%         end
%         end
%     end
    
    rfcall = [rfcall; rfc];
    pumall = [pumall; pum];
    pepall = [pepall; pep];
    
     
    
%     % For single, rotated mosaic            
%     [tok remain]=strtok(fNames{i},'_');
%     [tok remain]=strtok(remain,'_');
% 
%     refangle = str2double(tok);
%     
%     rfc = RFC{i}(RFC{i}~=-100);
%     pum = PUM{i}(PUM{i}~=-100);
%     pep = PEPPE{i}(PEPPE{i}~=-100);
%     
%     tmp = abs(rfc-refangle);
%     tmp(tmp>30) = 60-tmp( tmp>30 ); % If it's further than 30, subtract from 60 (phase unwrapping)
%     rfcabsdist{i} = tmp;
%     rfcabserr(i) = mean( tmp );
%     
%     tmp = abs(pum-refangle);
%     tmp(tmp>30) = 60-tmp( tmp>30 ); % If it's further than 30, subtract from 60 (phase unwrapping)
%     pumabsdist{i} = tmp;
%     pumabserr(i) = mean( tmp );
%     
%     tmp = abs(pep-refangle);
%     tmp(tmp>30) = 60-tmp( tmp>30 ); % If it's further than 30, subtract from 60 (phase unwrapping)
%     pepabsdist{i} = tmp;
%     pepabserr(i) = mean( tmp );

    end
end
% 
% rfcarray=[];
% pumarray=[];
% peparray=[];
% for i=1:length(rfcabsdist)
%     rfcarray = [rfcarray; rfcabsdist{i}];
%     pumarray = [pumarray; pumabsdist{i}];
%     peparray = [peparray; pepabsdist{i}];
% end
%     
% percentrfcright = (1- (sum(rfcarray(:)>3)/length(rfcarray(:))))*100
% percentpumright = (1- (sum(pumarray(:)>3)/length(pumarray(:))))*100
% percentpepright = (1- (sum(peparray(:)>3)/length(peparray(:))))*100
% 
% avgrfcerr = mean(rfcabserr)
% avgpumerr = mean(pumabserr)
% avgpeperr = mean(pepabserr)

% save 'SubMosaic_final';
% save 'PerfectRotatedMosaic_final'
% save 'AOSLO_Orientation';