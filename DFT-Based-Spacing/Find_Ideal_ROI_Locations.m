% Robert F Cooper 2018-11-15
%
% This script uses outputs from Montage_DFT_Analysis.m to find regions of
% interest.

clear;
close all force;

pName = '/local_data/Dropbox/Conference_Presentations/ARVO_2019/DFT_Methods/ROI_Test_Data';
fName = 'NC_11049_20160726_OD_confocal_Fouriest_Result.mat';


%[fName, pName]=

load(fullfile(pName,fName),'blendedim', 'blendederrim','threshold',...
                           'scaling', 'fovea_coords' );

nogozone = isnan(blendederrim);
blendederrim(isnan(blendederrim)) = 0;

blurerrim = imgaussfilt(blendederrim,64);

blurerrim(isnan(blurerrim))=1;


%%

x0 = [5.050, 2.195];

f = @(x)errfun(x,1-blurerrim, 16);
minmaxd = @(x)distfun(x,x0, 2, 256);

optim=optimoptions(@fmincon,'Display','iter','DiffMinChange',0.001,'ConstraintTolerance',1e-9);

[x, fval, exitflag]=fmincon(f, x0,[],[],[],[],[1 1], size(blurerrim),minmaxd,optim)


x0 = x0.*1000;
x = x.*1000;
figure; imagesc(1-blurerrim); colormap(jet(256)); axis image;
hold on;plot(x0(2),x0(1),'b*'); plot(x(2),x(1),'g*'); hold off;

function f=errfun(x, costim, roisize)
    x = x.*1000;

    halfroisize = roisize/2;
    roiranger = round( ((x(1)-halfroisize):(x(1)+halfroisize)) );
    roirangec = round( ((x(2)-halfroisize):(x(2)+halfroisize)) );

    f=mean2(costim(roiranger,roirangec));
end

function [c,ceq]=distfun(x, startpoint, mindist, maxdist)
    x = round(x.*1000);
    startpoint = round(startpoint.*1000);
    
    c = sqrt(sum((x-startpoint).^2))-maxdist; % Distance function

    
    ceq = double(sqrt(sum((x-startpoint).^2))==0);
end


