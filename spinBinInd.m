function [ spinBinInd ] = spinBinInd( spin_vec,M,calcSpinBin)
% Find the state index of the current state 'spin_vec' in the matrix M

[~,spinBinInd] = ismember(spin_vec',M,'rows');


end

