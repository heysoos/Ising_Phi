%% This code is used to calculate the KS test statistics between the distribution of the empirical data and the generalized data as a function of the temperature.

% KS test statistic calculates the maximum distance between two cumulative
% distributions.

clear all

load corr_all.mat
load fMRI_sub17.mat

n = 83;

for j = 1:10; % To calculate the test statistic for 10 realizations seperately
    corr_new = corr_DTI(:,:,:,j); % corr_DTI (83*83*252*10)
for i = 1:252;
    corr = corr_new(:,:,i);
    corr1 = reshape(corr,[n*n,1]);
    corr2 = reshape(Corr_FMRI,[n*n,1]);
    [~,~,ks2stat(j,i)] = kstest2(corr1,corr2); % Obtaining the test statistic for the 10 realizations as a function of temperature
end
end

% % figure
% 
% for k = 1:10;
%     figure
%     scatter(temp,ks2stat(k,:));
% end

% figure

% ks_mean = mean(ks2stat,1);
% scatter(temp,ks_mean)
