function [statistics] = radon_orientation_analysis(im, coords, um_per_pix)

imsize = size(im);

% Calculate the cell distances
mean_nn_dist=calc_icd(coords,um_per_pix);

micronboxsize = mean_nn_dist.*4.5;

pixelboxsize = ceil(micronboxsize./um_per_pix);

if any(pixelboxsize > imsize)    
    pixelboxsize = imsize(1);
end

pixelboxhalf = floor(pixelboxsize/2);

% center it on 0 so that we can rotate/translate to the center

anglecoords = [-3, 0;
                3, 0];

[X, Y]=meshgrid(0:(pixelboxhalf*2), 0:(pixelboxhalf*2));

X = X-pixelboxhalf;
Y = Y-pixelboxhalf;

mask = uint8(X.^2 + Y.^2 < pixelboxhalf.^2);

[V,C] = voronoin(coords,{'QJ'});

% Determine the number of sides for each coordinate location
for i=1:length(C)
   
    vertices=V(C{i},:);
    

    if (all(C{i}~=1)  && all( vertices(:,1)<=max(coords(:,1))) && all(vertices(:,2)<=max(coords(:,2)) ) ... % Second row are the Column limits
                      && all( vertices(:,1)>=min(coords(:,1))) && all(vertices(:,2)>=min(coords(:,2)) ) )

        numedges(i)=size(V(C{i},1),1);

    end
 
end

alist=-100*ones(size(coords,1),1);
expected_angle=-100*ones(size(coords,1),1);
for i=1:size(coords,1)
        
    x = round(coords(i,1));
    y = round(coords(i,2));
    
    if (y-pixelboxhalf) > 0           && (x-pixelboxhalf) > 0 && ...
       (y+pixelboxhalf) <= size(im,1) && (x+pixelboxhalf) <= size(im,2) ...
       && numedges(i)~=0
       
        roi = im( y-pixelboxhalf : y+pixelboxhalf, ...
                  x-pixelboxhalf : x+pixelboxhalf);

       
        roi = roi.*mask;
        outroi = double(roi.*mask);

        expected_angle(i) = (180*(numedges(i)-2))/numedges(i);
        range{i} = (90 - (expected_angle(i)/4)) : (90 + (expected_angle(i)/4));
        
        radonim = radon(roi,range{i})';

        clear numpeaks;
        
        % Determine the cutoffs for the radon xform by looking for the FWHM
        gausfilt = fspecial('gaussian',[1 5],.75);
        middle_orient = round( size(radonim,1)/2);
        
        middlerow = conv(radonim(middle_orient,:),gausfilt,'valid');
        
        xings = 1:length(middlerow);
        xings = xings(middlerow > (max(radonim(middle_orient,:))/2) );

        rrms=[];
        for r=1:size(radonim,1)
            radonrow = radonim(r,:);
            
            
            radonrow = conv(radonrow,gausfilt,'valid');
            radonrow = radonrow( (xings(1)):(xings(end)));
            dervline =  ( diff( radonrow,2) ); 
            
            rrms(r) = rms( dervline );

        end 
        
        rrms= rrms';

                  
        [rmsval angle] = max(rrms);
        [val worstangle] = min(rrms);
        
        theta = (angle - (1+(expected_angle(i)/4)) );

        alist(i) = theta;
        arms(i) = rmsval;
    end
end



patchsize = pixelboxsize;
patchwidths = 0:patchsize:imsize;

kk=1;
for ii=1:(length(patchwidths)-1)
beginpatchy = patchwidths(ii);
endpatchy   = ceil(patchwidths(ii+1));

    for jj=1:(length(patchwidths)-1)
        beginpatchx = patchwidths(jj);
        endpatchx   = ceil(patchwidths(jj+1));


        [patchcoords{kk} patches{kk}] = coordclipv2(coords, [beginpatchx endpatchx], [beginpatchy endpatchy]);
        kk=kk+1;
    end
end

for ii=1:length(patches)
    opatches{ii} = alist(patches{ii});
    bad_angles = opatches{ii} ~= -100;
    
    opatches{ii} = opatches{ii}(bad_angles);
    opatch_exp_angle{ii} = expected_angle(patches{ii});
    opatch_exp_angle{ii} = opatch_exp_angle{ii}(bad_angles);
end

variation=[];
weight=[];
for i=1:length(opatches) 
    current_patch = opatches{i};
    current_exp_angle = opatch_exp_angle{i};
    if length( current_patch ) > 4
        current_patch = current_patch(current_patch ~=-100);
        % Phase unwrap because we're doing std dev and need direct
        % distances
        [maxdiff maxind] = max(current_patch);

        alldists = pdist2(current_patch, current_patch(maxind));

        wrapping = abs( alldists ) > (current_exp_angle/4);
        if any( wrapping )
            current_patch( wrapping ) = current_patch( wrapping ) + current_exp_angle(wrapping)/2;
            current_patch( ~wrapping ) = current_patch( ~wrapping ) + alldists( ~wrapping ).*2;
        end
        
        variation = [variation var(current_patch)];
        weight = [weight length(current_patch)];
    end
end

top=0;
bottom=0;
for i=1:length(variation)
     top = top + ((weight(i)-1) * variation(i));
     bottom = bottom + weight(i);
end

statistics.Orientation_Pooled_Variance = top/bottom;
