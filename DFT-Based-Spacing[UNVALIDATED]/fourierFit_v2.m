function [maxnegdiff_ind ] = fourierFit_v2(fourierProfile)


%% Set up initial guess for fit parameters
doplots = false;

% Remove any nan and inf.
fourierProfile = fourierProfile(~isnan(fourierProfile));
fourierProfile = fourierProfile(~isinf(fourierProfile));
fourierProfile = fourierProfile-min(fourierProfile);

timeBase = 0:(length(fourierProfile)-1);

% Plot
if doplots
    thePlot = figure(10); clf; hold on
    set(gca,'FontName','Helvetica','FontSize',14);
    plot(fourierProfile,'k');
end


% Make initial guesses    
fitParams.scale1 = max(fourierProfile)*0.9-min(fourierProfile);
fitParams.decay1 = .05;
fitParams.offset1 = 0;
fitParams.shift = 0;


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

vlb = [0.01 0.001  -10  1];
vub = [15   15     10    length(fourierProfile)-2];

x = fmincon(@(x)FitModelErrorFunction(x,timeBase,fourierProfile,fitParams),x1,[],[],[],[],vlb,vub,[],options);

% Extract fit parameters
fitParams = XToParams(x,fitParams);

% Add final fit to plot
predictions = ComputeModelPreds(fitParams,timeBase);

% spacing = fitParams.shift;
if doplots
    figure(thePlot); hold on; plot(timeBase,predictions,'g','LineWidth',2); 
    hold off;drawnow;
end

% Find the second zero crossing (where the fit intersects with the curve)
residuals = predictions-fourierProfile;

residuals = medfilt1(residuals,7);

[pks,locs] = findpeaks(residuals*-1);
locs = locs(pks>0); % Find all local minima that are below 0 (the fit is underneath the data)

curheight = 1;
curind=1;

prevals = residuals(1:3);
for l=1:length(locs) % For each of the minima underneath the data,
        
    for i=locs(l):length(residuals)-1

        thisval = residuals(i);

        if all(prevals<0) && thisval>0 % Find the zero crossings
            curind=i;
            break;
        end
        prevals(1:2) = prevals(2:3);
        prevals(3) = thisval;
    end
    % If the zero crossing was preceded by a lower minima,
    % then take/keep that as our starting point for the next step.
    if residuals(locs(l))<=curheight
        curheight = residuals(locs(l));
        maxnegdiff_ind = locs(l);
    end
end

if doplots
    figure(11); clf;
    hold on; plot( maxnegdiff_ind, residuals(maxnegdiff_ind),'b*' );
end

% Trace back to where it is maximally different from our fit.
preval = residuals(maxnegdiff_ind-1)-residuals(maxnegdiff_ind);
for i=maxnegdiff_ind-1:-1:2
   
    thisval = residuals(i-1)-residuals(i);
    
    if round(preval, 5)<=0 && round(thisval,5)<=0 % It should only be decreasing or flat- if it isn't anymore and heads upward, kick out.
        maxnegdiff_ind=i; 
    elseif thisval>0.03
        if doplots
            figure(thePlot); hold on;
            plot( maxnegdiff_ind, fourierProfile(maxnegdiff_ind),'r*' );
            hold off;
        end
        break;
    end
    preval = thisval;
end



if doplots
    figure(11); plot( residuals );
    hold on; plot( maxnegdiff_ind, residuals(maxnegdiff_ind),'r*' );
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
    x = [params.scale1 params.decay1 params.offset1 params.shift];
end


% fitParams = XToParams(x,params)
%
% Convert search params and base structure to filled in structure.
function params = XToParams(x,params)
params.scale1 = x(1);
params.decay1 = x(2);
params.offset1 = x(3);
% params.scale2 = x(4);
% params.decay2 = x(5);
params.shift = x(4);
end

% preds =  ComputeModelPreds(params,t)
%
% Compute the predictions of the model
function fullExp = ComputeModelPreds(params,freqBase)

fullExp = params.offset1 + params.scale1*exp( -params.decay1 * (freqBase-params.shift) );

% bottomExpLoc = find(freqBase>params.shift);
% bottomExpTime = freqBase(bottomExpLoc);
% params.shift
% % The exponential must always line up with the other exponential function's
% % value!   
% maxmatch = fullExp(bottomExpTime(1))-params.scale2;
% 
% fullExp(bottomExpLoc) =  params.scale2*exp( -params.decay2 * (bottomExpTime-bottomExpTime(1)) );

end
