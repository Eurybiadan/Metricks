function [spacing, predictions, err, fitParams] = fourierFit(fourierProfile, prior, doplots)

if ~exist('doplots')
    doplots = false;
end


%% Set up initial guess for fit parameters

% Remove any nan and inf.
fourierProfile = fourierProfile(~isnan(fourierProfile));
fourierProfile = fourierProfile(~isinf(fourierProfile));
fourierProfile = fourierProfile-min(fourierProfile);
timeBase = 1:(length(fourierProfile));

fourierSampling =(1:length(fourierProfile))/(size(fourierProfile,2)*2);

%% Start plot
if doplots
    thePlot = figure(1); clf; hold on
    set(gca,'FontName','Helvetica','FontSize',14);
    plot(timeBase, fourierProfile,'k');
end

if isempty(prior)
    
    [fitParams.shift, firsterr] = fourierFit_v2(fourierProfile, doplots);
    % Make initial guesses
    fitParams.scale1 = 1;
    fitParams.decay1 = (fourierProfile(1)*.36) /...
                        (fitParams.shift-1);

    [maxval, maxind] = max(fourierProfile);
    if maxind ~= 1 % If the maximum value isn't the first index, 
                   % then ensure that the fit doesn't start touching the
                   % data
        maxval = maxval+1;
    end
                    
    fitParams.offset1 = maxval-fitParams.scale1;
    fitParams.scale2 =  fitParams.offset1*.3679;
    fitParams.decay2 = (fourierProfile(fitParams.shift)*.36) /...
                        (length(fourierProfile)-fitParams.shift);
        
else
    fitParams = prior;
end

% Add initial guess to the plot
predictions0 = ComputeModelPreds(fitParams,timeBase);
if doplots
    figure(thePlot); hold on; plot(timeBase,predictions0,'k','LineWidth',2); hold off;
end

%% Fit

% Set fmincon options
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','interior-point');

x1 = ParamsToX(fitParams);

vlb = [0.5 0.001 0.01 0.001  0.001  10];
vub = [5 0.5   15   15     0.5    length(fourierProfile)-2];

x = fmincon(@(x)FitModelErrorFunction(x,timeBase,fourierProfile,fitParams),x1,[],[],[],[],vlb,vub,[],options);

% Extract fit parameters
fitParams = XToParams(x,fitParams);

% Add final fit to plot
predictions = ComputeModelPreds(fitParams,timeBase);

if doplots
    figure(thePlot); hold on; plot(timeBase,predictions,'g','LineWidth',2);
    axis([0 150 0 5]);
end


residuals = fourierProfile-predictions;
spacing = ceil(fitParams.shift);

if doplots
    figure(2); clf; plot(residuals); hold on; plot(medfilt1(residuals,3));
    plot(spacing, residuals(spacing),'b*'); 
end

residuals = medfilt1(residuals,3);
preval = residuals(spacing-1)-residuals(spacing);

%%
minbound = 4;
maxbound = length(fourierProfile)-2;

for i=spacing-1:-1:minbound
   
    thisval = residuals(i-1)-residuals(i);
    
    if preval>=0 && thisval>=-0.01 % It should only be increasing or flat- if it isn't anymore and heads down, kick out.
        spacing=i; 

    elseif thisval<-0.01 && ((residuals(i-1)>0) || (residuals(i)>0))
        spacing=i;
        if doplots
            figure(thePlot); 
            plot(spacing, fourierProfile(spacing),'r*')
        end        
        break;
    end
    preval = thisval;
end

%% Determine Sharpness of the peak as an error measurment
lowfreqbound=spacing;
highfreqbound=spacing;

% f = fit([1:length(residuals)]',(fourierProfile-predictions)','smoothingspline','SmoothingParam', 0.3);
sharpresiduals = residuals; %f(1:length(residuals))';
% if doplots
%     figure(2); plot(sharpresiduals);
% end

%% Use a smoothed residual to find the bottoms of our peaks.
for i=(spacing-1):-1:minbound 
   
    thisval = sharpresiduals(i-1)-sharpresiduals(i);
    
    if thisval<=0.01 
        lowfreqbound=i; 

    elseif thisval>0.01
        lowfreqbound=i; 
        if doplots
            figure(2); hold on;
            plot(lowfreqbound, residuals(lowfreqbound),'g*')
        end
        break;
    end
    preval = thisval;
end
%%
for i=(spacing+1):1:maxbound
   
    thisval = sharpresiduals(i+1)-sharpresiduals(i);
    
    if thisval<=0.01 
        highfreqbound=i; 

    elseif thisval>0.01
        highfreqbound=i;
        if doplots
            figure(2); hold on;
            plot(highfreqbound, residuals(highfreqbound),'g*')
        end
        break;
    end
    preval = thisval;
end

maxamplitude = max(residuals(minbound:maxbound))-min(residuals(minbound:maxbound));

if lowfreqbound==(spacing-1) && highfreqbound~=spacing
    
    highheight = (residuals(spacing) - residuals(highfreqbound));
    highrun = fourierSampling(highfreqbound)-fourierSampling(spacing);

    heightdistinct = highheight./maxamplitude;
    
elseif highfreqbound==(spacing+1) && lowfreqbound~=spacing
    
    lowheight = (residuals(spacing) - residuals(lowfreqbound));
    lowrun = fourierSampling(spacing)-fourierSampling(lowfreqbound);

    heightdistinct = lowheight./maxamplitude;
    
elseif highfreqbound~=(spacing+1) && lowfreqbound~=(spacing-1)
    % Find the distinctness of our peak based on the average height of the two
    % sides of the triangle
    lowheight = residuals(spacing) - residuals(lowfreqbound);
    highheight = residuals(spacing) - residuals(highfreqbound);
    
    lowrun = fourierSampling(spacing)-fourierSampling(lowfreqbound);
    highrun = fourierSampling(highfreqbound)-fourierSampling(spacing);

    avgheight = (lowheight+highheight)/2;
%     avgrun = (lowrun+highrun)/2;

    heightdistinct = max([lowheight highheight])./maxamplitude;
else
    heightdistinct=0;
end


% Coefficient of determination
% SSres = sum(residuals.^2);
% SStot = sum( (fourierProfile - mean(fourierProfile)).^2 );
% n = length(fourierProfile);
% p = length(x)-1;

% err = 1 - ( (SSres./(n-p-1)) ./ (SStot./(n-1)) );

% err = sum(residuals(2:end).^2);

% spacing_ratio = (length(fourierProfile)./spacing);

err =  heightdistinct; %(err/firsterr); 

if doplots
    
    figure(2);
    hold on; plot(spacing, residuals(spacing),'r*');
    hold off;
    figure(1); title([' Quality: ' num2str(err) ]);
        hold off;
    drawnow;
%     pause;
end



end

% f = FitModelErrorFunction(x,timeBase,theResponse,fitParams)
%
% Search error function
function f = FitModelErrorFunction(x,timeBase,theResponse,fitParams)

% Extract parameters into meaningful structure
fitParams = XToParams(x,fitParams);

% Make predictions
preds = ComputeModelPreds(fitParams,timeBase);

% Compute fit error as RMSE
nPoints = length(theResponse);
theDiff2 = (theResponse-preds).^2;
f = 100*sqrt(sum(theDiff2)/nPoints);
% figure(333); hold on; plot(f,'.'); hold off;
end

% x = ParamsToX(params)
%
% Convert parameter structure to vector of parameters to search over
function x = ParamsToX(params)
    x = [params.scale1 params.decay1 params.offset1 params.scale2 params.decay2 params.shift];
end


% fitParams = XToParams(x,params)
%
% Convert search params and base structure to filled in structure.
function params = XToParams(x,params)
params.scale1 = x(1);
params.decay1 = x(2);
params.offset1 = x(3);
params.scale2 = x(4);
params.decay2 = x(5);
params.shift = x(6);
end

% preds =  ComputeModelPreds(params,t)
%
% Compute the predictions of the model
function fullExp = ComputeModelPreds(params,freqBase)

fullExp = params.offset1 + params.scale1*exp( -params.decay1 * freqBase );

bottomExpLoc = find(freqBase>params.shift);
bottomExpTime = freqBase(bottomExpLoc);

% The exponential must always line up with the other exponential function's
% value!   
maxmatch = fullExp(bottomExpTime(1))-params.scale2;

fullExp(bottomExpLoc) = maxmatch + params.scale2*exp( -params.decay2 * (bottomExpTime-bottomExpTime(1)) );

end
