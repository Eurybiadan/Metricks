function [spacing,orientation,relativeweight6th] = FFT2SpacingOrientationExctraction(Image,um_per_pix)
% Summary of the function goes here
%   Detailed explanation goes here

Npixel = size(Image,1);             % Ipothesis that image is a square matrix which Npixel x Npixel;
dim    = nextpow2(Npixel);          % Dimension of the image is power of 2; 
xmax   = 2^(dim-1)*um_per_pix;      % image half dimension in pixel (i.e. um_per_pix=1) or micron: i.e. image goes from -xmax to +xmax;
dx     = 2*xmax/Npixel;             % dimension of a pixel
x      = -xmax:dx:(xmax-dx);        

IMfft2 = fft2(Image);                         % FFT2 evaluation
deltaF = 1/(2*xmax); fmax = 1/(2*dx);         % units: cyc/micron or cyc/pixel (i.e um_per_pix)
freq   = -fmax:deltaF:(fmax - deltaF);

% IMfft2 = fft2(Image,2^(dim+1),2^(dim+1));       % FFT2 evaluation
% NNpixel = size(IMfft2,1);
% xxmax  = NNpixel/2;
% dxx    = 2*xxmax/NNpixel;
% deltaF = 1/(2*xxmax); fmax = 1/(2*dxx);           % units: cyc/micron or cyc/pixel (i.e um_per_pix)
% freq   = -fmax:deltaF:(fmax - deltaF);

Pow = abs(IMfft2).^2/Npixel^2;                     % Average of PSD for each frequency components; %max(max(abs(IMfft2).^2));% Normalized Power Spectrum
% Pow = abs(IMfft2).^2/NNpixel^2;                     % Average of PSD for each frequency components; %max(max(abs(IMfft2).^2));% Normalized Power Spectrum

[u,v] = meshgrid(freq,freq);



% figure(1)
% surfc(u,v,(fftshift(Pow))/max(max((fftshift(Pow)))),'FaceColor','interp','EdgeColor','none','FaceLighting','phong');%fftshift(Pow-Pow(1,1))
% axis tight, box on, grid on, axis square, view(2);%, colorbar;
% 
% imwrite( 255*(fftshift(Pow))/max(max(fftshift(Pow))), parula(256), 'Fourier_step.tif')



% Polar transformation: RHO & THETA
[theta,rho] = cart2pol(u,v);
Powfftshift = fftshift(Pow);

%Definition of Interpolation Power Spectra Function
Pow_Iterp = scatteredInterpolant(rho(:),theta(:),Powfftshift(:),'nearest');

% RHO   = 0: max(max(rho))/Npixel:max(max(rho)) - max(max(rho))/Npixel;
% rhoinit = 15;    %index for the first element of RHO to take into account: we discarge the DC value elements of Power Spectra
RHO   = 0.05: (max(max(rho))-0.05)/Npixel:max(max(rho)) - (max(max(rho))-0.05)/Npixel;
rhoinit = 1;    %index for the first element of RHO to take into account: we discarge the DC value elements of Power Spectra
THETA = -pi:2*pi/Npixel:pi-2*pi/Npixel;
% THETA = -pi:2*pi/NNpixel:pi-2*pi/NNpixel;

[qrho,qtheta] = meshgrid(RHO,THETA);
qPow = Pow_Iterp(qrho,qtheta);

% figure(2)
% surfc(qrho(:,rhoinit:end),qtheta(:,rhoinit:end),qPow(:,rhoinit:end)./max(max(qPow(:,rhoinit:end))),'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
% xlabel '\rho'; ylabel '\theta';
% axis tight, box on, grid on, axis square
% colorbar; 
% axis([0.05 0.5 min(THETA) max(THETA)])
% view(2)

% imwrite( 255*qPow./max(qPow(:)), parula(256), 'Pseudopolar_step.tif')

% Calculation of 1D FFT for only theta content of the Power Spectrum
for i = 1:size(qPow,2)
    
    qPowFT1D(:,i) = (fftshift(fft(abs(qPow(:,i)).^2)));
    
%     figure(3); imagesc(abs(qPowFT1D).^2);
end

% PowFT1D_MOD   = abs(qPowFT1D)'./max(max(abs(qPowFT1D))); % extraction of the modulo: look now we take the transpose of qPowFT1D
% PowFT1D_VOL   = sum(sum(abs(qPowFT1D')));                % pseudo-area under ALL PSD 
% PowFT1D_MOD   = abs(qPowFT1D')./PowFT1D_VOL;             % Modulo Normalization: force all PSDs to sum to one 

% PowFT1D_VOL   = sum(sum(abs(qPowFT1D(4:8,:)')));         % pseudo-area under 4th to 8 th Omega frequency PSD components 
PowFT1D_MOD   = abs(qPowFT1D').^2;%./PowFT1D_VOL;              


POWFT1D_PHASE = angle(qPowFT1D');                        % extraction of the phase

OMEGA = -Npixel/2:1:Npixel/2-1;                         

% figure(3)
% imagesc(OMEGA,RHO(rhoinit:end),PowFT1D_MOD(rhoinit:end,:)/max(max(PowFT1D_MOD(rhoinit:end,:))));
% set(gca, 'ydir', 'normal' );
% xlabel '\omega'; ylabel '\rho'; 
% axis tight, box on, grid on, axis square
% axis([1 10 0.1 0.3])
% title('Modulo');
% colorbar; 
% pause

% imwrite( 255*PowFT1D_MOD(rhoinit:end,:)'/max(max(PowFT1D_MOD(rhoinit:end,:))), parula(256), 'Modulus_step.tif')

% figure(4)
% imagesc(OMEGA,RHO,PowFT1D_MOD);
% set(gca, 'ydir', 'normal' )
% xlabel '\omega'; ylabel 'rho'; 
% title('Phase');
% colorbar; 
% axis([1 10 0.1 0.3])

% Post processing: extraction of the 6th frequency and its maximum value RHO (spacing) and
% its corresponding angle (orientation)

frequency_6th = find(OMEGA == 6);

distancefreq = PowFT1D_MOD(:,frequency_6th)';
% 
% frequency_4th = find(OMEGA == 4);
% frequency_5th = find(OMEGA == 5);
% frequency_7th = find(OMEGA == 7);
% frequency_8th = find(OMEGA == 8);

% weight6th = sum(PowFT1D_MOD(:,frequency_6th)');
% weight4th = sum(PowFT1D_MOD(:,frequency_4th)');
% weight5th = sum(PowFT1D_MOD(:,frequency_5th)');
% weight7th = sum(PowFT1D_MOD(:,frequency_7th)');
% weight8th = sum(PowFT1D_MOD(:,frequency_8th)');
% 
% relativeweight6th = weight6th./(weight4th+weight5th+weight6th...
%     +weight7th+weight8th);

% figure(5)
% plot(RHO(rhoinit:end),distancefreq(rhoinit:end)./max(distancefreq(rhoinit:end)),'r');
% xlabel '\rho '; ylabel 'maximum value normalizated [a.u]';
% grid on, axis square;
% saveas(gcf,'omega6plot.eps','epsc');

[weight6th,idx] = max(PowFT1D_MOD(rhoinit:end,frequency_6th));

% [weight4th,dummy] = max(PowFT1D_MOD(rhoinit:end,frequency_4th));
% [weight5th,dummy] = max(PowFT1D_MOD(rhoinit:end,frequency_5th));
% [weight7th,dummy] = max(PowFT1D_MOD(rhoinit:end,frequency_7th));
% [weight8th,dummy] = max(PowFT1D_MOD(rhoinit:end,frequency_8th));

% [weight4th,dummy] = max(PowFT1D_MOD(idx,frequency_4th));
% [weight5th,dummy] = max(PowFT1D_MOD(idx,frequency_5th));
% [weight7th,dummy] = max(PowFT1D_MOD(idx,frequency_7th));
% [weight8th,dummy] = max(PowFT1D_MOD(idx,frequency_8th));
% % 
% relativeweight6th = weight6th./(weight4th+weight5th+weight6th...
%     +weight7th+weight8th);

weight6th_all = sum(PowFT1D_MOD(rhoinit:end,frequency_6th)');

relativeweight6th = weight6th/weight6th_all; % value near 1 meaning that the cone arrangement is hexagonal, 
                                             % whereas value near 0 meaning poor hexagonal arrangment

exagon_size  = 1/RHO(rhoinit+idx-1);
exagon_phase = POWFT1D_PHASE(rhoinit+idx-1,frequency_6th);  %rad

spacing         = exagon_size;             % spacing of the cones
orientation     = -exagon_phase/6;         % orientation of the cones

% Computation of Spectral texture metrics - not included yet
% Srad = sum(qPow,1);
% Sang = sum(qPow,2);
% 
% S    = abs(fftshift(IMfft2)); 
% S    = mat2gray(log(1+S));
% 
% figure(5)
% imagesc(freq,freq,S);
% axis square, axis tight, box on; 
% grid on,   colorbar;
% title('Log Spectrum Normalizated');
% 
% figure(6)
% plot(RHO,Srad,'r-','Marker','.','MarkerSize',24,'LineWidth',6);
% xlabel '\rho '; ylabel 'Srad texture metric';
% grid on, axis square;
% 
% figure(7)
% plot(THETA*180/pi,Sang,'b-','Marker','.','MarkerSize',24,'LineWidth',6);
% xlabel '\theta '; ylabel 'Sang texture metric';
% grid on, axis square


end

