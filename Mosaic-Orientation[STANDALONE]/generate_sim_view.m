function [ sim_out, aligned_coords ] = generate_sim_view( cCoords )
% Robert Cooper 06-18-2012
%   This program generates the simulated view of the mosaic, generated from
%   rod and cone coordinates and masks.
%
%   Edited 04-19-2015
%   Set distribution to what we see in normal retina (ARVO 2014)
%
    scal = 10;
    
    [Y,X] = meshgrid(1:30*scal,1:30*scal);
    stddev_size=3.75*scal;
    cMask=exp(-(( X - ((30*scal)/2) )/stddev_size).^2-((Y- ((30*scal)/2) )/stddev_size).^2);


    cCoords = [cCoords(:,1)-min(cCoords(:,1)) cCoords(:,2)-min(cCoords(:,2))];
    cCoords = (scal*cCoords*[0 1;1 0])+1;
    

    cPad = 0; %size(cMask,2)*2;
    rPad = 0; %size(cMask,1)*2;

    
    sim_out = zeros( max(round(cCoords(:,1))) + round(rPad/2), round(max(cCoords(:,2))) + round(cPad/2) );

    shiftedcCoords = [round(cCoords(:,1)) + round(rPad/4), round(cCoords(:,2)) + round(cPad/4)];
    % Shift all coordinates to indicie space
    coordind = sub2ind(size(sim_out), shiftedcCoords(:,1), shiftedcCoords(:,2));
    
    % Generate the intensity range of the photoreceptors
    
    m=80; % Based off of target mean
    v=24.9^2; % Based off of known distributions (Cooper, ARVO 2014)
    
    mu = log((m^2)/sqrt(v+m^2));
    sigma = sqrt(log(v/(m^2)+1));
    
    int_values = lognrnd(mu,sigma,length(coordind),1);
%     int_values = 40;
    
    int_values(int_values>255) = 255;
    
    % Insert the random intensties
    sim_out(coordind) = uint8(int_values);
        
    % Perform the convolution to populate the simulation
    sim_out = imresize( conv2(sim_out,cMask,'same'), 1/scal) ;

%     minx = min(shiftedcCoords(:,1));
%     miny = min(shiftedcCoords(:,2));
%     shiftedcCoords = [shiftedcCoords(:,1)-minx shiftedcCoords(:,2)-miny]/scal;
%     shiftedcCoords = [shiftedcCoords(:,1)+minx shiftedcCoords(:,2)+miny];
    aligned_coords=((shiftedcCoords)*[0 1;1 0]-1)/scal;

end




