function [spacing, predictions, fitParams] = fourierFit_beta(fourierProfile, prior)

doplots = true;


%% Set up initial guess for fit parameters

% Remove any nan and inf.
fourierProfile = fourierProfile(~isnan(fourierProfile));
fourierProfile = fourierProfile(~isinf(fourierProfile));
fourierProfile = fourierProfile-min(fourierProfile);
timeBase = 0:(length(fourierProfile)-1);

fourierProfile = fourierProfile(2:end);
timeBase = timeBase(2:end);

%% Start plot
if doplots
    thePlot = figure(1); clf; hold on
    set(gca,'FontName','Helvetica','FontSize',14);
    plot(fourierProfile,'k');
end

if isempty(prior)
    
    fitParams.shift = fourierFit_v2(fourierProfile);
    shiftlb = fitParams.shift-25;
    shiftub = fitParams.shift+25;
    
    % Make initial guesses
    fitParams.scale1 = fourierProfile(2);
%     fitParams.decay1 = (fourierProfile(1)*.36) /...
%                         (fitParams.shift-1);
%     fitParams.exponent1 = exp(1);
%     fitParams.scale2 =  fourierProfile(fitParams.shift)*.3679;
%     fitParams.decay2 = (fourierProfile(fitParams.shift)*.36) /...
%                         (length(fourierProfile)-fitParams.shift);
%     fitParams.exponent2 = exp(1);
    
    fitParams.decay1 = -log(fourierProfile(round(fitParams.shift))./fitParams.scale1)./...
                        (fitParams.shift);

    fitParams.exponent1 = exp(1);
    fitParams.scale2 =  fourierProfile(round(fitParams.shift));
    fitParams.decay2 = -log(fourierProfile(140)./fourierProfile(round(fitParams.shift)))./...
                        (140);
    fitParams.exponent2 = exp(1);
        
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
options = optimset(options,'Diagnostics','off','Display','on','LargeScale','off',...
                   'Algorithm','sqp', 'MaxFunEvals', 10000);

x1 = ParamsToX(fitParams);

vlb = [0.5 0.001 1  1  0.001 1  shiftlb];
vub = [5   0.5   15 10 0.5   15 shiftub];

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
residuals = medfilt1(residuals,3);

preval = residuals(spacing-1)-residuals(spacing);
figure(2);
plot(spacing, residuals(spacing),'b*');

for i=spacing-1:-1:2
   
    thisval = residuals(i-1)-residuals(i);
    
    if preval>=0 && thisval>=0 % It should only be increasing or flat- if it isn't anymore and heads down, kick out.
        spacing=i; 

    elseif thisval<0.07
        if doplots
            figure(thePlot); 
            plot(spacing, fourierProfile(spacing),'r*')
        end
        break;
    end
    preval = thisval;
end

% [pks, locs]= findpeaks( residuals(1:ceil(spacing)) );
% 
% % If the last point is rising, then take add it to our list.
% if residuals(ceil(spacing))-residuals(ceil(spacing-1)) > 0
%     pks = [pks residuals(ceil(spacing))];
%     locs = ceil(spacing);
% end
% 
% if ~isempty(locs)
%     
%     pks = fliplr(pks);
%     locs = fliplr(locs);
%     if doplots
%         plot(locs(1), fourierProfile(locs(1)),'r*')
%     end
%     spacing = locs(1);
% end

if doplots
    hold off;drawnow;
    figure(2);hold on; plot(residuals); hold on; plot(spacing, residuals(spacing),'r*');
    hold off;
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
    x = [params.scale1 params.decay1 params.exponent1 params.scale2 params.decay2 params.exponent2 params.shift];
end


% fitParams = XToParams(x,params)
%
% Convert search params and base structure to filled in structure.
function params = XToParams(x,params)
params.scale1 = x(1);
params.decay1 = x(2);
params.exponent1 = x(3);
params.scale2 = x(4);
params.decay2 = x(5);
params.exponent2 = x(6);
params.shift = x(7);
end

% preds =  ComputeModelPreds(params,t)
%
% Compute the predictions of the model
function fullExp = ComputeModelPreds(params,freqBase)

fullExp = params.scale1*params.exponent1.^( -params.decay1 * freqBase );

bottomExpLoc = find(freqBase>=params.shift);
bottomExpTime = freqBase(bottomExpLoc);

% The exponential must always line up with the other exponential function's
% value!   
maxmatch = fullExp(bottomExpLoc(1))-params.scale2;

fullExp(bottomExpLoc) = maxmatch + params.scale2*params.exponent2.^( -params.decay2 * (bottomExpTime-bottomExpTime(1)) );

end
