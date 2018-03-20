%=========================================================================
%
% Monte Carlo Methods for the evaluation of the thermodynamics properties
% of the one-dimensional Ising model without taking the nearest-neighbour
% coupling
% ========================================================================
% method = 'glaub' or 'met'
% for glauber/metropolis transitions
%% Simulation Parameters
function  [ Ener,Mag, Spec_Heat, Sus, S, temp, Corr_DTI, J, time, fileSave, pet, Phi ] = Ising_Par_Phi_met_glaub(J, filename, calcPhi, method)

tic
N = length(J); % Dimension of the connectivity matrix

% temp = linspace(0.01,5,1000); %0.00001:0.001:3;
temp = logspace(-1,log10(4),200);
time_corr = 500;

thermalizeTime = 0.3; % How many percent of no_flip should you skip to let system thermalize?
no_flip = 100*N^2; % No of flip at each temperature
no_bin = 1; % # of bin

avgDenom = no_flip*(1-thermalizeTime); % Used for averaging the thermodynamic properties
lenTemp = length(temp);


saveFileBool = true;
fileSave = filename;

%% Connectivity Matrix

pet = ones(1,N);

% ------------------------------------------------------------------------
%% Matrices for Phi

if calcPhi == true
    calcTPM = true;
    calcSpinBin = true;
else
    calcTPM = false;
    calcSpinBin = false;
end

M = fliplr(((dec2bin(0:(2^N)-1)=='0') - 0)*2 - 1);
% Mtest = fliplr((dec2bin(0:(2^N)-1)) == '1');
dE_sys = M*J.*M.*2;

% ------------------------------------------------------------------------
% Generate Partition Indicies

PartIndex = 1:N;
% Generates all possible combinations of partitions from size 2:(MAX-1)
for i = 1:(size(PartIndex,2)-1)
    combinationsSet{i} = nchoosek(PartIndex,i);
    % NOTE: nchoosek function valid when length(PartIndex) < 15.
end
clear PartIndex

numPartSet = size(combinationsSet,2);

%% Montecarlo Simulation
TPM = zeros(2^N,N,lenTemp);
Phi = zeros(1,lenTemp);
spinBin = zeros(2^N,lenTemp);

Ener = zeros(1,lenTemp);
Mag = zeros(1,lenTemp);
Sus = zeros(1,lenTemp);
Spec_Heat = zeros(1,lenTemp);
S = zeros(N,time_corr,lenTemp);
Corr_DTI = zeros(N,N,lenTemp);

fprintf('Progress:\n');
fprintf(['\n' repmat('.',1,lenTemp) '\n\n']);
parfor tT = 1:lenTemp
    T = temp(tT);
    fprintf('\b|\n');
    
    spin_vec = zeros(1,N);
    randcoord = randperm(N);
    
    %%%IIT: TPM
    TPM(:,:,tT) = fpm_ising_met_g_sbn(dE_sys,M,J,T,calcTPM,method);
    phiMean = 0;
    spinBinMean = zeros(2^N,1);
   
    %%%
    
    enerMean = 0;
    enerSqrMean = 0;
    magMean = 0;
    magSqrMean = 0;
    susMean = 0;
    specHeatMean = 0;
    
    for b = 1 :no_bin
        spin_vec = (rand(N,1) > 0.5)*2 - 1;
                
        phiSum = 0;
        spinBinSum = zeros(2^N,1);
        
        enerSum = 0;
        enerSqrSum = 0;
        magSum = 0;
        magSqrSum = 0;
        nn = 0;
        
        for i = 1:no_flip
            % Run through the N points of the coordinate matrix once, and
            % when the counter is larger than N, it goes back to the first
            % point, and the coordinate is rerandomized before it is called
            % upon. This ensures the random but unique choice of index when
            % checking to flip.
            nn = nn + 1;
            if (nn > N)
                nn = nn-N;
                randcoord = randperm(N);
            end
            Flip = randcoord(nn);
            % Compute the change in energy
            dE=0;
            for j=1:N
                if (j~= Flip)
                    dE = dE + J(Flip,j)*spin_vec(j);
                end
            end
            dE=2*dE*spin_vec(Flip);
            
            if strcmp(method,'met')
                %%%%%%% METROPOLIS %%%%%%%%
                if (dE <= 0)
                    spin_vec(Flip) = - spin_vec(Flip);
                elseif (rand <= exp(-dE*pet(Flip)/T))
                    spin_vec(Flip) = - spin_vec(Flip);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            elseif strcmp(method,'glaub')
                
                %%%%%%%%% GLAUBER %%%%%%%%%
                if ( rand <= ( 1 / ( 1 + exp(dE/T) )) )
                    spin_vec(Flip) = - spin_vec(Flip);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            else
                error('Please choose an approprite method: glaub/met.')
            end
            
            if (i/no_flip) > thermalizeTime
                ener = 0;
                for ii = 1:N
                    for jj=ii+1:N
                        ener = ener - J(ii,jj)*spin_vec(ii)*spin_vec(jj);
                    end
                end
                
                %%%IIT: PHI!
                % phi = phi_ising_vir(numPartSet,combinationsSet,Mtest,TPM(:,:,tT),spin_vec,N,calcPhi);
                phi = 0;
                
                stateInd = state_ind(spin_vec,M,calcSpinBin);
                spinBinSum(stateInd) = spinBinSum(stateInd) + 1;
                
                phiSum = phiSum + phi;
                %%%
                
                enerSum = enerSum + ener;
                enerSqrSum = enerSqrSum + ener^2;
                magSum = magSum + abs(sum(spin_vec));
                magSqrSum = magSqrSum + sum(spin_vec)^2;
                
            end
            
        end
        
        phiMean(tT,b) = phiSum/avgDenom;
        spinBinMean(:,tT,b) = spinBinSum./avgDenom;
        
        enerMean(tT,b) = enerSum/avgDenom;
        magMean(tT,b) = magSum/avgDenom;
        enerSqrMean(tT,b) = enerSqrSum/avgDenom;
        magSqrMean(tT,b) = magSqrSum/avgDenom;
        specHeatMean(tT,b) = (enerSqrMean(tT,b)-(enerMean(tT,b))^2)/N/T^2;
        susMean(tT,b) = (magSqrMean(tT,b)-(magMean(tT,b))^2)/N/T;
        
    end
    Phi(tT) = mean(phiMean(tT,:)); %Divide by N? Doesn't make that much sense, but it might be a useful measure if we're interested in scaling effects...
    spinBin(:,tT) = mean(spinBinMean(:,tT,:),3);
    
    Ener(tT) = mean(enerMean(tT,:))/N;
    Mag(tT) = mean(magMean(tT,:))/N;
    Sus(tT)= mean(susMean(tT,:))/N;
    Spec_Heat(tT) = mean(specHeatMean(tT,:))/N;
    
    mm = 0;
    STemp = zeros(N,time_corr,lenTemp);
    for jjCorr = 1:time_corr
        for i = 1:N
            % Compute the change in energy
            mm = mm + 1;
            if (mm > N)
                mm = mm-N;
                randcoord = randperm(N);
            end
            Flip = randcoord(mm);
            dE = 0;
            for j=1:N
                if (j~=Flip)
                    dE = dE + J(Flip,j)*spin_vec(j);
                end
            end
            dE=2*dE*spin_vec(Flip);
            % Decide whether the change is possible
            if (dE <= 0)
                spin_vec(Flip) = - spin_vec(Flip);
            elseif (rand <= exp(-dE*pet(Flip)/T))
                spin_vec(Flip) = - spin_vec(Flip);
            end
            
        end
        STemp(:,jjCorr,tT)= spin_vec;  % temperature x timepoints x spinsites
    end
    S(:,:,tT) = STemp(:,:,tT);
    % Correlation matrix
    Corr_DTI(:,:,tT)=corrcoef(squeeze(S(:,:,tT))');
    
end

%%
clearvars -except Phi TPM spinBin Ener Mag Spec_Heat Sus S temp Corr_DTI J time fileDir fileSave saveFileBool 

if saveFileBool 
    save(fileSave)
end

toc
time = toc;
