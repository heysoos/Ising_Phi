%% This code is used to calculate the KS test statistics between the distribution of the empirical data and the generalized data as a function of the temperature.

% KS test statistic calculates the maximum distance between two cumulative
% distributions.

clear all

load Corr_DTI_6x6.mat
load Corr_FMRI_DMN_6x6.mat
% load DMN_6x6.mat

% Corr_FMRI = Corr_FMRI(DMN_6x6,DMN_6x6);
Corr_DTI(isnan(Corr_DTI)) = 0;
Corr_FMRI = Corr_FMRI2;
n = 6;

for i = 1:size(Corr_DTI,3);
    corr = Corr_DTI(:,:,i);
    corr1 = reshape(corr,[n*n,1]);
    corr2 = reshape(Corr_FMRI,[n*n,1]);
    [~,~,ks2stat(i)] = kstest2(corr1,corr2); % Obtaining the test statistic for the 10 realizations as a function of temperature
end


% for j = 1:10; % To calculate the test statistic for 10 realizations seperately
%     corr_new = corr_DTI(:,:,:,j); % corr_DTI (83*83*252*10)
% for i = 1:252;
%     corr = corr_new(:,:,i);
%     corr1 = reshape(corr,[n*n,1]);
%     corr2 = reshape(Corr_FMRI,[n*n,1]);
%     [~,~,ks2stat(j,i)] = kstest2(corr1,corr2); % Obtaining the test statistic for the 10 realizations as a function of temperature
% end
% end

% % figure
% 
% for k = 1:10;
%     figure
%     scatter(temp,ks2stat(k,:));
% end

% figure

% ks_mean = mean(ks2stat,1);
% scatter(temp,ks_mean)
