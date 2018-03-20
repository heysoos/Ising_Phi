function [ ks2stat, pDist2, corrInd ] = struct_corr_dist( dataStruct, fSrch, srchStr, fCorrDist, Corr_DTI )
%   struct_avg
%
%   Input: struct - struct to look at
%
%   fSrch - the field whose values we are searching/filtering through
%
%   srchStr - string which the field fSrch will be selecting for
% 
%   Output: ks2stat and pDist2 are the distances between the input Corr_DTI
%   (the simulated correlation matrix) with respect to fCorrDist (the
%   empirical correlation matrix)
%


for i = 1:numel(dataStruct) % use loop since dimensions for this string are not equal
    
    fSrchVal2{i,:} =  dataStruct(i).(fSrch);
    
end
corrInd = ismember(fSrchVal2,srchStr);
tempStruct = dataStruct(corrInd);

% tic
% [ks2stat, pDist2] = arrayfun(@(x) ks_test_dist(Corr_DTI, x.(fCorrDist)), tempStruct, 'UniformOutput', false);
% toc
tic
for i = 1:length(tempStruct)
   [ks2stat{i}, pDist2{i}] = ks_test_dist(Corr_DTI, tempStruct(i).(fCorrDist)); 
end
display(['Correlation distance calculations finished in ',num2str(toc), ' seconds'])
end

