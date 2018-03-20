%% N = 5 , Full, Random,  No Phi, Metropolis & Glauber %%%%%%%
tic
calcPhi = true;

for i = 1:200
    
    sparsity = 0;
    J = make_J(5,sparsity);
    J = J./max(max(J));
    
    % METROPOLIS
    method = 'met';
    
    dir = ['Simulations/Ising_random/N5_motif_full/', method,'/'];
    name = ['Ising_random_',method,'_',num2str(i),'.mat'];
    filename = [dir, name];
    
    Ising_Par_Phi_met_glaub(J,filename, calcPhi, method);
    
    % GLAUBER
    
    method = 'glaub';
    
    dir = ['Simulations/Ising_random/N5_motif_full/', method,'/'];
    name = ['Ising_random_',method,'_',num2str(i),'.mat'];
    filename = [dir, name];
    
    Ising_Par_Phi_met_glaub(J,filename, calcPhi, method);
    
end
tFinal = toc/60;

fprintf('Total computation time: %2.2f minutes', tFinal)
%% N = 5, Networks %%%%%%%%%%%%%%%%%%%%%
clc; clear; 

Networks = {'Aud', 'DMN', 'ECN_L', 'ECN_R', 'Salience', 'Sensorimotor', 'VISL', 'VISM', 'VISO'};
load('/home/heysoos/MATLAB-Drive/Phi/reduced_networks/indices_84.mat') % used to extract J

numIters = 1;

dirN = 'Simulations/Ising_Networks/';
calcPhi = true;

% Jtemp = load('MJ_all.mat');
% Jtemp = Jtemp.MJ_all_mn;

Jtemp = load('/home/heysoos/MATLAB-Drive/Phi/connectome_data/mean_struct_corr.mat');
Jtemp = Jtemp.meanJ_prob;

for iNet = 1:length(Networks)
    
    Ind = indices.(Networks{iNet});
    
    dist = distances(graph(1./Jtemp));
    J = 1./dist(Ind,Ind);
    J(isinf(J)) = 0;
    J = J./max(max(J));
    
    for i = 1:numIters
        
        % METROPOLIS
        method = 'met';
        
        wd = [dirN, Networks{iNet},'/',method,'/'];
        if ~exist(wd,'dir')
            mkdir(wd)
        end
        
        name = ['Ising_', method,'_',num2str(i),'.mat'];
        filename = [wd, name];
                
        Ising_Par_Phi_met_glaub(J,filename, calcPhi, method);
        
        % GLAUBER
        method = 'glaub';
        
        wd = [dirN, Networks{iNet},'/',method,'/'];
        if ~exist(wd,'dir')
            mkdir(wd)
        end
        
        name = ['Ising_', method,'_',num2str(i),'.mat'];
        filename = [wd, name];
                
        Ising_Par_Phi_met_glaub(J,filename, calcPhi, method);
        
    end
end

% LEGACY CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% N = 5 , Full, Random,  Phi %%%%%%%%%%%%%%%%%%%%%
% tic
% dir = 'Simulations/Ising_Phi_random/N5_motif_full_newPhi_vir/';
% calcPhi = true;
% 
% for i = 1:200
%     
%     name = ['Ising_Phi_random_',num2str(i),'.mat'];
%     filename = [dir, name];
%     
%     sparsity = 0;
%     J = make_J(5,sparsity);
%     J = J./max(max(J));
%     
%     % Ising_Par_Phi(J,filename, calcPhi);
%     Ising_Par_Phi_vir(J,filename, calcPhi);
%     
% end
% tFinal = toc/60;
% 
% fprintf('Total computation time: %2.2f minutes', tFinal)
%% N = 5, Networks %%%%%%%%%%%%%%%%%%%%%
% Networks = {'Aud', 'DMN', 'ECN_L', 'ECN_R', 'Salience', 'Sensorimotor', 'VISL', 'VISM', 'VISO'};
% load('/home/heysoos/Dropbox/Masters Thesis/Codes/Phi/reduced_networks/indices.mat') % used to extract J
% 
% dirN = 'Simulations/Ising_Phi_Networks/';
% calcPhi = true;
% 
% Jtemp = load('MJ_all.mat');
% Jtemp = Jtemp.MJ_all_mn;
% 
% for iNet = 1:length(Networks)
%     
%     for i = 1:200
%         
%         wd = [dirN, Networks{iNet},'/'];
%         if ~exist(wd,'dir')
%             mkdir(wd)
%         end
%         name = ['Ising_Phi',num2str(i),'.mat'];
%         filename = [wd, name];
%                 
%         Ind = indices.(Networks{iNet});
%         
%         % J = J./mean(mean(J));
%         
%         dist = distances(graph(1./Jtemp));
%         J = 1./dist(Ind,Ind);
%         J(isinf(J)) = 0;
%         J = J./max(max(J));
%         
%         Ising_Par_Phi(J,filename, calcPhi);
%         
%     end
% end

%% N = 7, random %%%%%%%%%%%%%%%%%%%%%
% dir = 'Simulations/Ising_Phi_random/N7';
% 
% 
% for i = 1:200
%     
%     name = ['Ising_Phi_random_',num2str(i),'.mat'];
%     filename = [dir, name];
%     
%     J = make_J(7,rand);
%     Ising_Par_Phi(J,filename, calcPhi);
%     
% end