

clear;
close all;

% [fname, pname] = uigetfile(pwd);

fname = 'Test_Group_results.csv';
pname = pwd;

data = dlmread( fullfile(pname, fname) );

eccent_raw = data(:,1);
icd_spac = data(:,2);
dft_spac = data(:,3);

eccent_range = [150 250;
                300 425;
                425 550;
                800 1000;
                1200 1700];
            
simulations = 50000;

count = zeros(size(eccent_range,1),1);
eccent = cell(size(eccent_range,1),1);
icd_spac_eccent = cell(size(eccent_range,1),1);
dft_spac_eccent = cell(size(eccent_range,1),1);
            
for i=1:size(eccent_range,1)
    
    over = eccent_raw >= eccent_range(i,1);
    under = eccent_raw <= eccent_range(i,2);
    
    whats_the = over & under;
    
    eccent{i} = eccent_raw(whats_the);
    icd_spac_eccent{i} = icd_spac(whats_the);
    dft_spac_eccent{i} = dft_spac(whats_the);
    count(i) = sum(whats_the);
end

rng('shuffle')
rando_grid = zeros(simulations,length(count));

for j=1:size(eccent_range,1)

    rando_grid(:,j) = randi([1 count(j)], simulations, 1);
    
end

fit_eccent = zeros(size(rando_grid));
fit_icd = zeros(size(rando_grid));
fit_dft = zeros(size(rando_grid));



for m=1:simulations    
    for i=1:size(fit_eccent,2)
        fit_eccent(m,i) = eccent{i}( rando_grid(m,i) ) ;
        fit_icd(m,i) = icd_spac_eccent{i}( rando_grid(m,i) ) ;
        fit_dft(m,i) =  dft_spac_eccent{i}( rando_grid(m,i) ) ;
    end
end
icd_power = zeros(3,simulations);
icd_residuals = zeros(size(eccent_range,1),simulations);
dft_power = zeros(3,simulations);

parfor m=1:simulations
    
    [the_fit, gof, fit_info] = fit(fit_eccent(m,:)', fit_icd(m,:)','power2');
    icd_power(:,m) = [the_fit.a; the_fit.b; the_fit.c];
    icd_residuals(:,m) = fit_info.residuals;
    
%     the_fit = fit(fit_eccent(m,:)', fit_dft(m,:)','power2');
%     dft_power(:,m) = [the_fit.a; the_fit.b; the_fit.c];

%     if the_fit.a > 100
%     figure(2); plot(the_fit)
%     hold on; 
%     plot(fit_eccent(m,:)', fit_dft(m,:)')
%     end
end

%% Plot everything

mean_icd_fit = mean(icd_power,2)
icd_confidence = 1.96*std(icd_power,[],2);

% mean_dft_fit = mean(dft_line,2);
% mean_dft_fit(1) = 10.^mean_dft_fit(1);
% dft_confidence = 1.96*std(dft_line,[],2);
% dft_confidence(1) = 10.^dft_confidence(1);

x = (1:max(eccent_raw)+100);

figure(1); histogram(icd_power(1,icd_power(1,:)>-20 & icd_power(1,:)<20) , 20)
figure(2); histogram(icd_power(2,icd_power(2,:)>-20 & icd_power(2,:)<20) , 20)
figure(3); histogram(icd_power(3,icd_power(3,:)>-20 & icd_power(3,:)<20) , 20)

%plot((eccent_raw), (icd_spac),'b.'); hold on;
% plot(x, .0228*x.^0.7454 + 1.3471,'b');
% plot(x, (mean_icd_fit(1)+icd_confidence(1))* x.^(mean_icd_fit(2) + icd_confidence(2)) ,'b');

% plot((eccent_raw), (dft_spac),'r.');
% plot(x, mean_dft_fit(1)* x.^mean_dft_fit(2),'r');