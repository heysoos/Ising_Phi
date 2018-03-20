function [ phi ] = phi_ising(numPartSet,combinationsSet,Mtest,TPM,spin_vec,N,calcBool)

%  TO DO: In IIT 3.0, the purview and the mechanism are two seperate
%  concepts. In this version of the code, they are the same. The partitions
%  looked at in this code seperate both the mechanism that constrains our
%  calculation (the state of the partition) as well as the purview (which
%  is the set of elements that defines the space of the cause repertoire).

% This is an important distinction because, apparently because (check out the
% supplementary material of IIT 3.0 on the topic of virtual elements) one
% needs to consider the correlations between an element outside of the
% partition and elements within the partition. For example, if element C is
% outside of the partition and has causal connections with both A and B,
% then marginalizing over C = 1,0 would lead to correlations in the
% probabilities for mechanism AB.

% To account for this, one can calculate the marginal probabilities for
% each mechanism seperately and then multiply the cause-repertoires
% together at the very end. See the function phi_ising_vir for the new
% version of phi_ising that deals with this.

if ~calcBool
    phi = 0;
    return
end
%           -------------- PHI ------------
pCount = 0;
currentState = bi2de(logical(spin_vec' == -1));

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
        
        
        %%%%%%%%%%%%%%% TPM REDUCTION ALGORITHM %%%%%%%%%%%%%%
        % For each row in the whole system's cause repertoire, generates
        % indices representing to which partition in the cause repertoire
        % they fall under. REPRESENTS HOW THE ROWS OF THE TPM SHOULD BE
        % AVERAGED.
        
        % TO DO: this is unnecessary, just count in decimal to begin with,
        % why use Mtest?
        indP1 = bi2de(Mtest(:,part1));
        indP2 = bi2de(Mtest(:,part2));
        
        % Logical array representing all whole system states that are equal
        % to the partitioned current state. REPRESENTS COLUMNS OF THE TPM
        % The Columns of the TPM that correspond to the partitioned current
        % state
        
        indTPV1 = bi2de(Mtest(:,part1)) == bi2de(spin_vec(part1)' == -1);
        indTPV2 = bi2de(Mtest(:,part2)) == bi2de(spin_vec(part2)' == -1);
        
        TPM1 = mean(TPM(:,indTPV1),2);
        TPM2 = mean(TPM(:,indTPV2),2);
        
        % Transition probabily vector (TPV): Cause repertoire
        % of first and second partition
        
        TPV1 = zeros(1,2^length(part1));
        TPV2 = zeros(1,2^length(part2));
        TPV = zeros(1,2^N);
        
        % TO DO: Consider pre-calculating the logical indices for this loop
        % and remove the loop altogether....or at least mean outside the
        % loop.
        
        
        for iTPV1 = 1:2^length(part1)
            TPV1(iTPV1) = mean(TPM1((indP1 == (iTPV1 - 1)),:));
        end
        TPV1 = TPV1./(sum(TPV1));
        
        for iTPV2 = 1:2^length(part2)
            TPV2(iTPV2) = mean(TPM2((indP2 == (iTPV2 - 1)),:));
        end
        TPV2 = TPV2./(sum(TPV2));
        
        % TO DO: Consider using an outerproduct method instead of this
        % loop.
        
        for iFinal = 1:2^N
            TPV(iFinal) = TPV1(indP1(iFinal)+1) * TPV2(indP2(iFinal)+1);
        end
        
        
        
        % Kullback-Leibler Divergance
        %----------------------------------------------------
        % % TO DO: IIT 3.0 uses EMD (Earth Mover's Distance). Might solve
        % Phi > N % problem when N = 2? Check this out.
        %-----------------------------------------------------
        
        PTPM = TPM(:,currentState+1)./sum(TPM(:,currentState+1));
        
        H2 = PTPM.*log2(PTPM./TPV');
        H2(isnan(H2)) = 0;
        H2 = sum(H2);
%         effInfo2(pCount) = H2;
        H2 = H2/(min([length(part1) length(part2)]));
        effInfoNorm2(pCount) = H2;
        
        
        % TO DO: This normalization doesn't exist in IIT 3.0.
        
    end
end
phi = min(effInfoNorm2);

end