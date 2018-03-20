function [ FPM ] = fpm_ising_met_g_sbn( dE_sys,M,J,T,calcBool,method)
N = length(J);

if ~calcBool
    %THERE HAS TO BE A BETTER WAY THAN THIS OMG IM SO EMBARASSED
    FPM = zeros(2^N,N);
    return
end

if strcmp(method,'glaub') % GLAUBER TRANSITIONS
    FPM = 1./(1 + exp(dE_sys/T) ); 
    
elseif strcmp(method,'met') % METROPOLIS TRANSITIONS
    detFlip = (dE_sys <= 0);
    FPM = dE_sys; FPM(detFlip) = 1; FPM(~detFlip) = exp(-FPM(~detFlip)/T);
    
    
    for iTPM = 1:2^N
        for jTPM = 1:2^N
            
            lgclP = logical(M(iTPM,:) - M(jTPM,:));
            if sum(lgclP) > 1
                TPM(iTPM,jTPM) = 0;
            elseif sum(lgclP) == 1
                FPMtemp = FPM(iTPM,:);
                TPM(iTPM,jTPM) = FPMtemp(lgclP)/N;
                
            else
                FPMtemp = FPM(iTPM,:);
                FPMtemp(~lgclP) = 1 - FPMtemp(~lgclP);
                % TPM(iTPM,jTPM) = prod(FPMtemp);
                TPM(iTPM,jTPM) = mean(FPMtemp);
            end
            if sum(sum(isnan(TPM))) > 0
                
                pause()
                display('sum(sum(isnan(TPM))) > 0! TPM calculation paused...')
            end
            
        end
    end
else
    error('Please choose an approprite method: glaub/met.')
    
end

