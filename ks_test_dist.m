function [ks2stat, pKS, pDist] = ks_test_dist(Corr_DTI,Corr_FMRI)
% This code is used to calculate the KS test statistics between the
% distribution of the empirical data and the generalized data as a function
% of the temperature.

% KS test statistic calculates the maximum distance between two cumulative
% distributions.

% TO DO: 
% - fix how NaNs are dealt with using a custom distance function 


Corr_DTI(isnan(Corr_DTI)) = 1; % we can't just do this!

% n = size(Corr_DTI,2);
Corr_FMRI(1:length(Corr_FMRI)+1:end) = 0;

numSamples = size(Corr_DTI,3);

ks2stat = zeros(1,numSamples);
pKS = zeros(1,numSamples);
pDist = zeros(1,numSamples);

for iSample = 1:numSamples
    corr = Corr_DTI(:,:,iSample);
    
%     corr1 = reshape(corr,[n*n,1]);
%     corr2 = reshape(Corr_FMRI,[n*n,1]);

    corr(1:length(corr)+1:end) = 0;
    corr1 = squareform(corr);
    corr2 = squareform(Corr_FMRI);
    
    %%%%%% FISHER TRANSFORM (13-12-2017) %%%%%%%
    
    corr1 = atanh(corr1);
    corr2 = atanh(corr2);
    try 
        [~,pKS(iSample),ks2stat(iSample)] = kstest2(corr1,corr2);
    catch
        pKS(iSample) = NaN;
        ks2stat(iSample) = NaN;
    end
    
    if nargout > 1
        % euclidean distance
        % pDist(iSample) = pdist([corr1; corr2]);
        pDist(iSample) = pdist([corr1; corr2],@naneucdist); % ignore NaNs 
    end
end


% for j = 1:10; % To calculate the test statistic for 10 realizations seperately
%     corr_new = corr_DTI(:,:,:,j); % corr_DTI (83*83*252*10)
% for i = 1:252;
%     corr = corr_new(:,:,i);
%     corr1 = reshape(corr,[n*n,1]);
%     corr2 = reshape(Corr_FMRI,[n*n,1]);
%     [~,~,ks2stat(j,i)] = kstest2(corr1,corr2);% Obtaining the test statistic for the 10 realizations as a function of temperature
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


% took this from the matlab pdist docs page
function D2 = naneucdist(XI,XJ)  
%NANEUCDIST Euclidean distance ignoring coordinates with NaNs
n = size(XI,2);
sqdx = (XI-XJ).^2;
nstar = sum(~isnan(sqdx),2); % Number of pairs that do not contain NaNs
nstar(nstar == 0) = NaN; % To return NaN if all pairs include NaNs
D2squared = nansum(sqdx,2).*n./nstar; % Correction for missing coordinates
D2 = sqrt(D2squared);
