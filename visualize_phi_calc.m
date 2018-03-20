% clearvars 
%% params
J = make_J(6,0.5);
temp = linspace(0.1,5,1000);
%%
N = length(J);
%% Matrices for Phi

calcPhi = true;
calcTPM = true;
calcSpinBin = true;

M = fliplr(((dec2bin(0:(2^N)-1)=='0') - 0)*2 - 1);
Mtest = fliplr((dec2bin(0:(2^N)-1)) == '1');
dE_sys = M*J.*M.*2;

TPM = zeros(2^N);
PF = zeros(length(M),1);

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

%% calculate TPM


% for tT = 1:length(temp)
%     T = temp(tT);
%     T
%     TPM_full(:,:,tT) = tpm_ising(dE_sys,M,J,T,calcTPM);
% end
% 
% for i = 1:size(TPM,3)
%     imagesc(TPM_full(:,:,i))
%     txt1 = sprintf('T = %1.2f', temp(i));
%     xlabel(txt1)
%     axis square
%     pause(0.01)
% end
%% Define a random state and a random temperature

spin_vec = (rand(1,N)*2 - 1)';
spin_vec(spin_vec>0) = 1;
spin_vec(spin_vec<=0) = -1;

spin_vec = [1 1 -1 -1 1 1]';

tT = 200;

TPM = TPM_full(:,:,tT);

%% Phi

%           -------------- PHI ------------
pCount = 0;
currentState = bi2de(logical(spin_vec' == -1));

% Generates the simplified TPM (transition probability matrix).
% This matrix determines the nature of a transition from 3
% groups for the Ising model: Deterministic flips towards
% current state, deterministic flips away from current state,
% and probabilistic flips in either direction. This matrix is
% used for calculating the transition probability matrix in the
% transP_Par2 function.

for i_part = 1:ceil(numPartSet/2)
    partLength = length(combinationsSet{i_part});
    if i_part == ceil(numPartSet/2) && mod(numPartSet,2) ~= 0
        jLoop = partLength/2;
    else
        jLoop = partLength;
    end
    
    % Generate Phi for this bi-partition pair
    for j_part = 1:jLoop
        pCount = pCount + 1;
        part1 = combinationsSet{i_part}(j_part,:);
        part2 = combinationsSet{numPartSet - (i_part - 1)}(partLength - (j_part - 1),:);
        
        TPV1 = zeros(1,2^length(part1));
        TPV2 = zeros(1,2^length(part2));
        
        %%%%%%%%%%%%%%% TPM REDUCTION ALGORITHM %%%%%%%%%%%%%%
        % For each row in the whole system's cause repertoire,
        % generates indices representing to which partition in
        % the cause repertoire they fall under. REPRESENTS HOW
        % THE ROWS OF THE TPM SHOULD BE AVERAGED.
        indP1 = bi2de(Mtest(:,part1));
        indP2 = bi2de(Mtest(:,part2));
        
        % Logical array representing all whole system states
        % that are equal to the partitioned current state.
        % REPRESENTS COLUMNS OF THE TPM
        % The Columns of the TPM that correspond to the
        % partitioned current state
        indTPV1 = bi2de(Mtest(:,part1)) == bi2de(spin_vec(part1)' == -1);
        indTPV2 = bi2de(Mtest(:,part2)) == bi2de(spin_vec(part2)' == -1);
        
        TPM1 = mean(TPM(:,indTPV1),2);
        TPM2 = mean(TPM(:,indTPV2),2);
        
        % Transition probabily vector (TPV): Cause repertoire
        % of first and second partition
        for iTPV1 = 1:2^length(part1)
            TPV1(iTPV1) = mean(TPM1((indP1 == (iTPV1 - 1)),:));
        end
        TPV1 = TPV1./(sum(TPV1));
        for iTPV2 = 1:2^length(part2)
            TPV2(iTPV2) = mean(TPM2((indP2 == (iTPV2 - 1)),:));
        end
        TPV2 = TPV2./(sum(TPV2));
        
        for iFinal = 1:2^N
            TPV(iFinal) = TPV1(indP1(iFinal)+1) * TPV2(indP2(iFinal)+1);
        end
        
        
        
        % Kullback-Leibler Divergance
        %----------------------------------------------------
        % % TO DO: IIT 3.0 uses EMD (Earth Mover's Distance).
        % Might solve Phi > N % problem when N = 2? Check this
        % out.
        %-----------------------------------------------------
        
        PTPM = TPM(:,currentState+1)./sum(TPM(:,currentState+1));
        
        H2 = PTPM.*log2(PTPM./TPV');
        H2(isnan(H2)) = 0;
        H2 = sum(H2);
        effInfo2(pCount) = H2;
        H2 = H2/(min([length(part1) length(part2)]));
        effInfoNorm2(pCount) = H2;
        
        
        % TO DO: This normalization doesn't exist in IIT 3.0.
        
    end
end
phi = min(effInfoNorm2);