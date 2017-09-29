% close all;
clear;

pathname = pwd;
[fNames] = read_folder_contents(pathname,'csv');

% Load the scaling file
[tmp, lutData] = load_scaling_file('/AOSLO_Repo1/Personnel/Rob_Cooper/Orientation_Data/ROIs/LUT.csv');
% [tmp, lutData] = load_scaling_file('M:\Documents\Grad school\Cell_Orientation_2015\LUT.csv');
% [tmp, lutData] = load_scaling_file('/remote_project_folders/Rob Cooper/Cell_Orientation_2015/LUT.csv');

numnames = length(fNames);
for n=1:numnames
    
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

%     jitter = .45*randn(size(coords));
%     coords = coords + jitter;
    
    if ~exist('im','var')
        [im, shiftcoords] = generate_sim_view( coords ); 
        im = uint8(im);

        coords = shiftcoords;
%         imagesc(im); colormap gray; axis image;
    end
    
    RFC_coord_ONLY_radon;
    RFC{n} = alist;
        
    rfc_pum_mosaic_orientation_v3;
    PUM{n} = six_sided_deg;

    Main_FFT2ConePacking_verJ2015;
    PEPPE{n} = peppedegs;
    
    
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
    
    toanalyze(n) = percentsixsided(n) > 65;

    %% Draw the angles on the image for each algorithm
    % Parse out which angles are common to all algorithms
%% For Agreement
    matches = logical((alist~=-100) & (six_sided_deg~=-100) & (peppedegs~=-100));
    
    rfc = (alist(matches));
    rfcmatchcoords = coords(matches,:);
    pum = (six_sided_deg(matches));
    pummatchcoords = coords(matches,:);
    pep = (peppedegs(matches));
    pepmatchcoords = coords(matches,:);
    unmatchcoords  = coords(~matches,:);
    
%% For match to perfect
%     rfc = (alist(alist~=-100));
%     rfcmatchcoords = coords(alist~=-100,:);
%     pum = (six_sided_deg(six_sided_deg~=-100));
%     pummatchcoords = coords(six_sided_deg~=-100,:);
%     pep = (peppedegs(peppedegs~=-100));
%     pepmatchcoords = coords(peppedegs~=-100,:);
    
%% For agreement
    [unwraprfc,unwrappum,unwrappep, dists]= unwrap_phase(rfc,pum,pep);
    
    % rfcpum rfcpep peppum
    dists = abs(dists);
    allunwrap = [unwraprfc,unwrappum,unwrappep];
    algagree = all(dists <= 6,2);
    someagree = (dists <= 6);
    twoagree = sum(someagree,2)==2;
    
    [mindist mindistind]=min(dists(twoagree,:),[], 2);
    
    jj=1;
    for ii=1:length(twoagree)
        if twoagree(ii)
            someagree(ii,:) = 0;
            someagree(ii,mindistind(jj)) = 1;
            jj=jj+1;
        end
    end
    
%     agreeorient = mean(allunwrap,2);
    
    agreepercent(n) = sum(algagree)./ length(algagree);
    
    anglecoords = [-3, 0;
                    3, 0];

                
%% For match to perfect         
%     [tok remain]=strtok(fNames{n},'_');
%     [tok remain]=strtok(remain,'_');
% 
%     refangle = -str2double(tok);
%     
%     rfcoffang = abs( rfc-(refangle) );
%     rfcoffang = rfcoffang<=3;
%     if any(~rfcoffang)
%        disp( [fNames{n} ' has a GREATER THAN 3 cell using RADON!'])
%     end
%     
%     pumoffang = abs( pum-(refangle) );
%     pumoffang = pumoffang<=3;
%     if any(~pumoffang)
%        disp( [fNames{n} ' has a GREATER THAN 3 cell using PUM!'])
%     end
%     
%     pepoffang = abs( pep-(refangle) );
%     pepoffang = pepoffang<=3;
%     if any(~pepoffang)
%        disp( [fNames{n} ' has a GREATER THAN 3 cell using FOURIER!'])
%     end
%                 

%% General display
    figure(1); imagesc(im); colormap gray; axis image; hold on;
%     plot(unmatchcoords(:,1), unmatchcoords(:,2), 'wo');
    
    for ii=1: length(rfc)

        x = (rfcmatchcoords(ii,1));
        y = (rfcmatchcoords(ii,2));
        
        
        if all( someagree(ii,:) )
            agreeorient = mean([unwraprfc(ii), unwrappum(ii), unwrappep(ii)]);
            agreetform = affine2d([cosd(agreeorient) -sind(agreeorient) 0; sind(agreeorient) cosd(agreeorient) 0; x y 1]);
            [agreerotx agreeroty]=transformPointsForward(agreetform,anglecoords(:,1),anglecoords(:,2));
        elseif someagree(ii,1)
            agreeorient = mean([unwraprfc(ii), unwrappum(ii)]);
            agreetform = affine2d([cosd(agreeorient) -sind(agreeorient) 0; sind(agreeorient) cosd(agreeorient) 0; x y 1]);
            [agreerotx agreeroty]=transformPointsForward(agreetform,anglecoords(:,1),anglecoords(:,2));
        elseif someagree(ii,2)           
            agreeorient = mean([unwraprfc(ii), unwrappep(ii)]);
            agreetform = affine2d([cosd(agreeorient) -sind(agreeorient) 0; sind(agreeorient) cosd(agreeorient) 0; x y 1]);
            [agreerotx agreeroty]=transformPointsForward(agreetform,anglecoords(:,1),anglecoords(:,2));
        end
        
        % Draw the orientations of each mosaic on the image
        rfctform = affine2d([cosd(rfc(ii)) -sind(rfc(ii)) 0; sind(rfc(ii)) cosd(rfc(ii)) 0; x y 1]);
        
        [rfcrotx rfcroty]=transformPointsForward(rfctform,anglecoords(:,1),anglecoords(:,2));
        
        
        figure(1); hold on;
        
        plot(x, y,'w.');        

%         if ii ==44
%            disp('44'); 
%         end
        
        if all( someagree(ii,:) )
           plot(agreerotx,agreeroty,'y');  
        elseif any(someagree(ii,1:2))
           plot(agreerotx,agreeroty,'y'); 
        else
           plot(rfcrotx,rfcroty,'r'); 
        end
%         if rfcoffang(ii)
%             plot(rfcrotx,rfcroty,'c');
%         else
%             plot(rfcrotx,rfcroty,'r');            
%         end        
        
        title(strrep(fname(1:end-4),'_',' '));
        hold off;
        
    end
    
%     figure(2); imagesc(im); colormap gray; axis image; hold on;
    for ii=1: length(pum)

        x = (pummatchcoords(ii,1));
        y = (pummatchcoords(ii,2));
        
        if someagree(ii,3)           
            agreeorient = mean([unwrappum(ii), unwrappep(ii)]);
            agreetform = affine2d([cosd(agreeorient) -sind(agreeorient) 0; sind(agreeorient) cosd(agreeorient) 0; x y 1]);
            [agreerotx agreeroty]=transformPointsForward(agreetform,anglecoords(:,1),anglecoords(:,2));
        end
        
        % Draw the orientations of each mosaic on the image
        pumtform = affine2d([cosd(pum(ii)) -sind(pum(ii)) 0; sind(pum(ii)) cosd(pum(ii)) 0; x y 1]);        

        [pumrotx pumroty]=transformPointsForward(pumtform,anglecoords(:,1),anglecoords(:,2));        

        figure(1); hold on;

%         plot(x, y,'w.');

        if all( someagree(ii,:) )
            %nop
        elseif someagree(ii,3)
            plot(agreerotx,agreeroty,'y');
        elseif someagree(ii,1)
            % nop
        else
            plot(pumrotx,pumroty,'g'); 
        end
        
%         if pumoffang(ii)
%             plot(pumrotx,pumroty,'g'); 
%         else
%             plot(pumrotx,pumroty,'r');
%         end
        
        title(strrep(fname(1:end-4),'_',' '));
        hold off;
        
    end
    
%     figure(3); imagesc(im); colormap gray; axis image; hold on;
    for ii=1: length(pep)

        x = (pepmatchcoords(ii,1));
        y = (pepmatchcoords(ii,2));
        
        % Draw the orientations of each mosaic on the image
        peptform = affine2d([cosd(pep(ii)) -sind(pep(ii)) 0; sind(pep(ii)) cosd(pep(ii)) 0; x y 1]);

        [peprotx peproty]=transformPointsForward(peptform,anglecoords(:,1),anglecoords(:,2));
        
        figure(1); hold on;

%         plot(x, y,'w.');

        if all( someagree(ii,:) )
            %nop
        elseif someagree(ii,3) || someagree(ii,2)
            %nop
        else
            plot(peprotx,peproty,'b');
        end
%         if pepoffang(ii)
%             plot(peprotx,peproty,'b'); 
%         else
%             plot(peprotx,peproty,'r'); 
%         end
        
        title(strrep(fname(1:end-4),'_',' '));
        hold off;
        
    end
    
%     out = input(['Write this? (Percent agreement: ' num2str(agreepercent(n)*100) ' ) ']);
%     out = input('Write this? ');
%     if out == 1        
%         figure(1);
%         saveas(gcf, [fname(1:end-4) '_radon.eps'],'epsc');
%         figure(2);
%         saveas(gcf, [fname(1:end-4) '_pum.eps'],'epsc');
%         figure(3);
%         saveas(gcf, [fname(1:end-4) '_fourier.eps'],'epsc');
%         figure(1);
%         saveas(gcf, [fname(1:end-4) '_agreement.eps'],'epsc');
%         imwrite(im, [fname(1:end-4) '.tif']);
%     end
    
    %%
    pause(0.01);
    clear im;
   
end

%%
rfcall = [];
pumall = [];
pepall = [];
meansndiffs=[];
% Analyze their differences.
for i=1:length(fNames)

    i
 
    % For comparing submosaics and real mosaics
    % Find all cases where they match
    matches = logical((RFC{i}~=-100) & (PUM{i}~=-100) & (PEPPE{i}~=-100));
    rfc = (RFC{i}(matches));
    pum = (PUM{i}(matches));
    pep = (PEPPE{i}(matches));
    
    diffrfcpum = (rfc-pum);
    diffrfcpep = (rfc-pep);
    diffpeppum = (pep-pum);
    
    diffs = [diffrfcpum diffrfcpep diffpeppum];
    phasewrapped = abs(diffs) > 30;
    
%     Determine our biggest difference- this will be our anchor
    [vals inds]=max( abs( diffs ),[],2);
    
    % Shift the difference to the anchor
    for j=1:length(inds)
        
        if phasewrapped(j, inds(j))
        switch inds(j)
            case 1 % between rfc and pum
                if sign(diffrfcpum(j)) == 1 % rfc is higher, and the anchor
                    pum(j) = pum(j) +60;
                    
                    % Flip pep over the
                    % anchor as well to maintain relative distances
                    if abs(diffrfcpep(j)) > 30
                        pep(j) = pep(j)+60;
                    else                                                
                        pep(j) = diffrfcpep(j).*2+pep(j);
                    end
                    
                else % pum higher, and is the anchor- bring rfc to it
                    rfc(j) = rfc(j) +60;
                    
                    % Flip pep over the
                    % anchor as well to maintain relative distances
                    if abs(diffpeppum(j)) > 30
                        pep(j) = pep(j)+60;
                    else
                        pep(j) = -diffpeppum(j).*2+pep(j);
                    end
                    
                end

            case 2 % between rfc and pep
                if sign(diffrfcpep(j)) == 1 % rfc is the anchor
                    pep(j) = pep(j) +60;
                                     
                    % Flip pep over the
                    % anchor as well to maintain relative distances
                    if abs(diffrfcpum(j)) > 30
                        pum(j) = pum(j)+60;
                    else                      
                        pum(j) = diffrfcpum(j).*2+pum(j);                        
                    end                   
                    
                else % pep is the anchor
                    rfc(j) = rfc(j) +60;
                    
                    % Flip pep over the
                    % anchor as well to maintain relative distances
                    if abs(diffpeppum(j)) > 30
                        pum(j) = pum(j)+60;
                    else
                        pum(j) = diffpeppum(j).*2+pum(j);
                    end 
                    
                end
            case 3 % between pep and pum
                if sign(diffpeppum(j)) == 1 % pep is the anchor
                    pum(j) = pum(j) +60;
                                       
                    % Flip pep over the
                    % anchor as well to maintain relative distances
                    if abs(diffrfcpep(j)) > 30
                        rfc(j) = rfc(j)+60;
                    else
                        rfc(j) = -diffrfcpep(j).*2+rfc(j);
                    end                   
                    
                else % pum is the anchor
                    pep(j) = pep(j) +60;
                    
                    % Flip pep over the
                    % anchor as well to maintain relative distances
                    if abs(diffrfcpum(j)) > 30
                        rfc(j) = rfc(j)+60;
                    else
                        rfc(j) = -diffrfcpum(j).*2+rfc(j);
                    end 
                    
                end
        end
        end
    end
%     
%     diffrfcpum = (rfc-pum);
%     diffrfcpep = (rfc-pep);
%     diffpeppum = (pep-pum);
%     
%     meansndiffs = [meansndiffs; [(rfc+pum)/2 diffrfcpum (rfc+pep)/2 diffrfcpep (pep+pum)/2 diffpeppum]];
    
%     rfcall = [rfcall; rfc];
%     pumall = [pumall; pum];
%     pepall = [pepall; pep];
    
     
    
    % For single, rotated mosaic            
    [tok remain]=strtok(fNames{i},'_');
    [tok remain]=strtok(remain,'_');

    refangle = -str2double(tok);
    
    rfc = RFC{i}(RFC{i}~=-100);
    pum = PUM{i}(PUM{i}~=-100);
    pep = PEPPE{i}(PEPPE{i}~=-100);
    
    tmp = abs(rfc-refangle);
    tmp(tmp>30) = 60-tmp( tmp>30 ); % If it's further than 30, subtract from 60 (phase unwrapping)
    rfcabsdist{i} = tmp;
    rfcabserr(i) = mean( tmp );
    
    tmp = abs(pum-refangle);
    tmp(tmp>30) = 60-tmp( tmp>30 ); % If it's further than 30, subtract from 60 (phase unwrapping)
    pumabsdist{i} = tmp;
    pumabserr(i) = mean( tmp );
    
    tmp = abs(pep-refangle);
    tmp(tmp>30) = 60-tmp( tmp>30 ); % If it's further than 30, subtract from 60 (phase unwrapping)
    pepabsdist{i} = tmp;
    pepabserr(i) = mean( tmp );

 
end

rfcarray=[];
pumarray=[];
peparray=[];
for i=1:length(rfcabsdist)
    rfcarray = [rfcarray; rfcabsdist{i}];
    pumarray = [pumarray; pumabsdist{i}];
    peparray = [peparray; pepabsdist{i}];
end
    
percentrfcright = (1- (sum(rfcarray(:)>3)/length(rfcarray(:))))*100
percentpumright = (1- (sum(pumarray(:)>3)/length(pumarray(:))))*100
percentpepright = (1- (sum(peparray(:)>3)/length(peparray(:))))*100

avgrfcerr = mean(rfcabserr)
avgpumerr = mean(pumabserr)
avgpeperr = mean(pepabserr)

stdrfcerr = std(rfcabserr)
stdpumerr = std(pumabserr)
stdpeperr = std(pepabserr)

% save 'SubMosaic_final';
% save 'PerfectRotatedMosaic_JITTERLESS_final'
% save 'AOSLO_Orientation';