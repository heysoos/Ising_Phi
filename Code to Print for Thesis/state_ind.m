function [ spinBinInd ] = state_ind( spin_vec,M,calcSpinBin)
% Find the state index of the current state 'spin_vec' in the matrix M
if ~calcSpinBin
    spinBinInd = 1;
    return
end

[~,spinBinInd] = ismember(spin_vec',M,'rows');
% currentState = bi2de(logical(spin_vec' == -1)) + 1; % SLOWER APPARENTLY!


end

