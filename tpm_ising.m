function [ TPM ] = tpm_ising( dE_sys,M,J,T,calcBool)
N = length(J); 

if ~calcBool
    %THERE HAS TO BE A BETTER WAY THAN THIS OMG IM SO EMBARASSED
    TPM = zeros(2^N,2^N);
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%% TPM ALGORITHM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This generates the Flip Probablity Matrix (FPM) for the system. It is
% written in the LOLI (Low-Order bits correspond to Low-Index nodes)
% convention as described in the IIT 3.0 python documentation
% (https://pythonhosted.org/pyphi/examples/2014paper.html). This matrix
% is then used to generate the Transition Probability Matrix (TPM) for
% each temperature.


detFlip = (dE_sys <= 0);
FPM = dE_sys; FPM(detFlip) = 1; FPM(~detFlip) = exp(-FPM(~detFlip)/T);
% FPM = 1-M.*tanh(M*J'/T); FPM = FPM./(max(max(FPM))); % GLAUBER TRANSITIONS


for iTPM = 1:2^N
    for jTPM = 1:2^N
        
        lgclP = logical(M(iTPM,:) - M(jTPM,:));
        FPMtemp = FPM(iTPM,:);
        FPMtemp(~lgclP) = 1 - FPMtemp(~lgclP);
        TPM(iTPM,jTPM) = prod(FPMtemp);
        if sum(sum(isnan(TPM))) > 0
            
            pause()
            display('sum(sum(isnan(TPM))) > 0! TPM calculation paused...')
        end
        
    end
end

end

