%% Set default bg figure color to white
set(0,'defaultfigurecolor',[1 1 1])

fontname = 'Helvetica';
set(0,'defaultaxesfontname',fontname);
set(0,'defaulttextfontname',fontname);

fontsize = 12;
set(0,'defaultaxesfontsize',fontsize);
set(0,'defaulttextfontsize',fontsize);


addpath('../') % add the Phi folder to path
addpath(genpath('../cbrewer'))
addpath(genpath('../npy-matlab-master'))
addpath(genpath('../kakearney-boundedline-pkg-32f2a1f'))
addpath(genpath('../altmany-export_fig-5be2ca4'))

%% Load network data (phi/no phi)

homePhi = '../'; % Should lead to main 'Phi' directory
% load('/home/heysoos/Dropbox/Masters Thesis/Codes/Phi/reduced_networks/indices.mat')
load([homePhi, 'reduced_networks/indices_84.mat'])
simFolder = {'Aud', 'DMN', 'ECN_L', 'ECN_R', 'Salience', 'Sensorimotor', 'VISL', 'VISM', 'VISO'};
methods = {'met'; 'glaub'}; %  Metropolis and Glauber Transitions

iSub_met = 1;
iSub_glaub = 1;
for iSimFolder = simFolder    
    
    for method = 1:length(methods) % 2 transition methods
        
        wdirPhi = ['Simulations/pyPhi/Ising_Networks/',iSimFolder{1},'/'];
        wdir = ['Simulations/Ising_Networks/Current_Sims/',iSimFolder{1},'/'];
        
        mStr = char(methods(method));
        
        wdir = [wdir, mStr, '/'];
        wdirPhi = [wdirPhi, mStr,'/'];
        
%         switch mStr
%             case 'met'
%                 
%                 wdir = [wdir, 'met/'];
%                 wdirPhi = [wdirPhi, 'met/'];
%                 
%             case 'glaub'
%                 
%                 wdir = [wdir, 'glaub/'];
%                 wdirPhi = [wdirPhi, 'glaub/'];
%                 
%             otherwise
%                 error('You fucked up!')
%         end
        
        dirStruct = dir(wdir);
        dirStructPhi = dir(wdirPhi);
        
        numFiles = sum(~cat(1,dirStruct(:).isdir));
        numFilesPhi = sum(~cat(1,dirStructPhi(:).isdir));
        for i = 1:numFiles
            
            namePhi = ['Ising_',mStr,'_', num2str(i),'_meanPhi.npz']; % *************
            name = ['Ising_',mStr, '_', num2str(i),'.mat'];
            
            
            filename = [wdir, name];
            filenamePhi = [wdirPhi, namePhi];
            
            tempLoad = load(filename,'J','temp','Sus','Corr_DTI','Phi','Mag','Ener','Spec_Heat','S');
            
            if numFilesPhi == 1
                unzip(filenamePhi,wdirPhi)
            end
            
            Phi = readNPY([wdirPhi, 'phiSum.npy']);
            PhiSus = readNPY([wdirPhi, 'phiSus.npy']);
            T_Phi = readNPY([wdirPhi, 'T2.npy']);
            
%             dirStructPhi = dir(wdirPhi); % this updates the dirStruct   
            
            
            switch mStr
                case 'met'
                    data_NET_met(iSub_met).N = length(tempLoad.J);
                    data_NET_met(iSub_met).Network = iSimFolder{1};
                    data_NET_met(iSub_met).J = tempLoad.J;
                    data_NET_met(iSub_met).S = tempLoad.S;
                    data_NET_met(iSub_met).temp = tempLoad.temp;
                    data_NET_met(iSub_met).Mag = tempLoad.Mag;
                    data_NET_met(iSub_met).Ener = tempLoad.Ener;
                    data_NET_met(iSub_met).Sus = tempLoad.Sus;
                    data_NET_met(iSub_met).Spec_Heat = tempLoad.Spec_Heat;
                    data_NET_met(iSub_met).Corr_DTI = tempLoad.Corr_DTI;
                    data_NET_met(iSub_met).Phi = Phi';
                    data_NET_met(iSub_met).PhiSus = PhiSus';
                    data_NET_met(iSub_met).tempPhi = T_Phi';
                    
                    iSub_met = iSub_met + 1;
                    
                case 'glaub'
                    data_NET_glaub(iSub_glaub).N = length(tempLoad.J);
                    data_NET_glaub(iSub_glaub).Network = iSimFolder{1};
                    data_NET_glaub(iSub_glaub).J = tempLoad.J;
                    data_NET_glaub(iSub_glaub).S = tempLoad.S;
                    data_NET_glaub(iSub_glaub).temp = tempLoad.temp;
                    data_NET_glaub(iSub_glaub).Mag = tempLoad.Mag;
                    data_NET_glaub(iSub_glaub).Spec_Heat = tempLoad.Spec_Heat;
                    data_NET_glaub(iSub_glaub).Ener = tempLoad.Ener;
                    data_NET_glaub(iSub_glaub).Sus = tempLoad.Sus;
                    data_NET_glaub(iSub_glaub).Corr_DTI = tempLoad.Corr_DTI;
                    data_NET_glaub(iSub_glaub).Phi = Phi;
                    data_NET_glaub(iSub_glaub).PhiSus = PhiSus;
                    data_NET_glaub(iSub_glaub).tempPhi = T_Phi';
                    
                    iSub_glaub = iSub_glaub + 1;
            end
            
            
            
            
        end
        
    end
end

%% Load random data (phi/no phi)

homePhi = '../';
% load('/home/heysoos/Dropbox/Masters Thesis/Codes/Phi/reduced_networks/indices.mat')
load([homePhi, 'reduced_networks/indices_84.mat'])
methods = {'met'; 'glaub'}; %  Metropolis and Glauber Transitions

iSub_met = 1;
iSub_glaub = 1;

for method = 1:length(methods) % 2 transition methods
    
    wdirPhi = ['Simulations/pyPhi/Ising_random/N5_motif_full/'];
    wdir = ['Simulations/Ising_random/N5_motif_full/'];
    
    mStr = char(methods(method));
    
    wdir = [wdir, mStr, '/'];
    wdirPhi = [wdirPhi, mStr,'/'];
        
    dirStruct = dir(wdir);
    dirStructPhi = dir([wdirPhi,'*.npz']);
    
    tempExDir = [wdirPhi, 'temp/']; %unzips npz file here temporarily, deletes folder afterwards
    if ~exist(tempExDir,'dir')
        mkdir(tempExDir)
    end    
    
    numFiles = sum(~cat(1,dirStruct(:).isdir));
    numFilesPhi = sum(~cat(1,dirStructPhi(:).isdir));
    
    for i = 1:numFilesPhi
        
%         namePhi = ['Ising_random_',mStr,'_', num2str(i),'_meanPhi.npz']; % *************
        namePhi = dirStructPhi(i).name;
        simNumber = strsplit(namePhi,'_');
        simNumber = simNumber{4}; % get the file number label to match .npy and .mat files
        name = ['Ising_random_',mStr, '_', num2str(simNumber),'.mat'];
        
        
        filename = [wdir, name];
        filenamePhi = [wdirPhi, namePhi];
        
        tempLoad = load(filename,'J','temp','Sus','Corr_DTI','Phi','Mag','temp','Ener','Spec_Heat','S');
        

        unzip(filenamePhi,tempExDir)

        
        Phi = readNPY([tempExDir, 'phiSum.npy']);
        PhiSus = readNPY([tempExDir, 'phiSus.npy']);
        T_Phi = readNPY([tempExDir, 'T2.npy']);
        
        delete([tempExDir, '*.npy'])
        
        %             dirStructPhi = dir(wdirPhi); % this updates the dirStruct
        
        
        switch mStr
            case 'met'
                data_RAND_met(iSub_met).N = length(tempLoad.J);
                data_RAND_met(iSub_met).J = tempLoad.J;
                data_RAND_met(iSub_met).S = tempLoad.S;
                data_RAND_met(iSub_met).temp = tempLoad.temp;
                data_RAND_met(iSub_met).Mag = tempLoad.Mag;
                data_RAND_met(iSub_met).Spec_Heat = tempLoad.Spec_Heat;
                data_RAND_met(iSub_met).Ener = tempLoad.Ener;
                data_RAND_met(iSub_met).Sus = tempLoad.Sus;
                data_RAND_met(iSub_met).Corr_DTI = tempLoad.Corr_DTI;
                data_RAND_met(iSub_met).Phi = Phi';
                data_RAND_met(iSub_met).PhiSus = PhiSus';
                data_RAND_met(iSub_met).temp = tempLoad.temp;
                data_RAND_met(iSub_met).tempPhi = T_Phi';
                
                iSub_met = iSub_met + 1;
                
            case 'glaub'
                data_RAND_glaub(iSub_glaub).N = length(tempLoad.J);
                data_RAND_glaub(iSub_glaub).J = tempLoad.J;
                data_RAND_glaub(iSub_glaub).S = tempLoad.S;
                data_RAND_glaub(iSub_glaub).temp = tempLoad.temp;
                data_RAND_glaub(iSub_glaub).Mag = tempLoad.Mag;
                data_RAND_glaub(iSub_glaub).Spec_Heat = tempLoad.Spec_Heat;
                data_RAND_glaub(iSub_glaub).Ener = tempLoad.Ener;
                data_RAND_glaub(iSub_glaub).Sus = tempLoad.Sus;
                data_RAND_glaub(iSub_glaub).Corr_DTI = tempLoad.Corr_DTI;
                data_RAND_glaub(iSub_glaub).Phi = Phi;
                data_RAND_glaub(iSub_glaub).PhiSus = PhiSus;
                data_RAND_glaub(iSub_glaub).tempPhi = T_Phi';
                
                iSub_glaub = iSub_glaub + 1;
        end
        
        
        
        
    end
    
end

%%
clearvars -except data_NET_met data_NET_glaub data_RAND_met data_RAND_glaub
%% Playing with the new data (NET)
% Glaub and Met visualized together for each network

data_TEMP_met = data_NET_met;
data_TEMP_glaub = data_NET_glaub;

for iNet = 1:length(data_TEMP_met)
    
%     T = data_TEMP_met(iNet).temp;
    T = data_TEMP_met(iNet).temp;
    
    subplot(2,2,1)
    
    SusMet = data_TEMP_met(iNet).Sus;     
    
    plot(T,SusMet,'.','MarkerSize',11)    
    axis square
    title('Metropolis')
    xlabel('Temperature')
    ylabel('\chi')
    
    subplot(2,2,2)
    SusGlaub = data_TEMP_glaub(iNet).Sus;
    plot(T,SusGlaub,'.','MarkerSize',11)
    axis square
    title('Glauber')
    xlabel('Temperature')
    ylabel('\chi')
    
    subplot(2,2,3)
    
    Phi = data_TEMP_met(iNet).Phi;
    plot(T,Phi,'.','MarkerSize',11)
    axis square
    title('Metropolis')
    xlabel('Temperature')
    ylabel('\Phi')
    
    subplot(2,2,4)
    Phi = data_TEMP_glaub(iNet).Phi;
    plot(T,Phi,'.','MarkerSize',11)
    axis square
    title('Glauber')
    xlabel('Temperature')
    ylabel('\Phi')
    
    Network = data_TEMP_met(iNet).Network;
    suptitle(Network)
    pause()
end
%% Playing with the new data (RANDOM)
% Glaub/Met visualized seperately 

data_TEMP_met = data_RAND_met;
data_TEMP_glaub = data_RAND_glaub;

for iNet = 1:length(data_TEMP_met)
    
    subplot(1,2,1)
    
    T = data_TEMP_met(iNet).temp;
    SusMet = data_TEMP_met(iNet).Sus;     
    
    plot(T,SusMet,'.','MarkerSize',11)    
    axis square
    xlabel('Temperature')
    ylabel('\chi')
    
    subplot(1,2,2)
    
    T = data_TEMP_met(iNet).tempPhi;
    Phi = data_TEMP_met(iNet).Phi;
    
    plot(T,Phi,'.','MarkerSize',11)
    axis square
    xlabel('Temperature')
    ylabel('\Phi')
    
    suptitle(['Metropolis, sim #:', num2str(iNet)])
    
    pause()
    
end
    
for iNet = 1:length(data_TEMP_glaub)
    
    subplot(1,2,1)
    
    T = data_TEMP_glaub(iNet).temp;
    SusGlaub = data_TEMP_glaub(iNet).Sus;
    
    plot(T,SusGlaub,'.','MarkerSize',11)
    axis square
    xlabel('Temperature')
    ylabel('\chi')    
    
    subplot(1,2,2)
    
    T = data_TEMP_glaub(iNet).tempPhi;
    Phi = data_TEMP_glaub(iNet).Phi;
    
    plot(T,Phi,'.','MarkerSize',11)
    axis square
    xlabel('Temperature')
    ylabel('\Phi')
    
    suptitle(['Glauber, sim #:', num2str(iNet)])
    
    pause()
end

%% Same plot as above but overlapped Sus/Phi



for iNet = 1:9
    
    T = data_NET_met(iNet).temp;   
    
    SusMet = data_NET_met(iNet).Sus;
    SusMet = SusMet./max(max(SusMet));
    
    SusGlaub = data_NET_glaub(iNet).Sus;
    SusGlaub = SusGlaub./max(max(SusGlaub));
    
    PhiMet = data_NET_met(iNet).Phi;
    PhiMet = PhiMet./max(max(PhiMet));
    PhiGlaub = data_NET_glaub(iNet).Phi;
    PhiGlaub = PhiGlaub./max(max(PhiGlaub));
    
    PhiSusMet = data_NET_met(iNet).PhiSus;
    PhiSusMet = PhiSusMet./max(max(PhiSusMet));
    PhiSusGlaub = data_NET_glaub(iNet).PhiSus;
    PhiSusGlaub = PhiSusGlaub./max(max(PhiSusGlaub));
    
    figure(1)
    subplot(2,2,1)
    
%     plot(T,SusMet,'.',T,PhiMet,'.',T,PhiSusMet,'.','MarkerSize',11)
    plot(T,SusMet,'.',T,PhiMet,'.','MarkerSize',11)
    legend({'\chi','\Phi'})
    axis square
    title('Metropolis')
    xlabel('Temperature')
    ylabel('\chi')
    
    subplot(2,2,2)
    
%     plot(T,SusGlaub,'.',T,PhiGlaub,'.',T,PhiSusGlaub,'.','MarkerSize',11)
    plot(T,SusGlaub,'.',T,PhiGlaub,'.','MarkerSize',11)
    legend({'\chi','\Phi'})
    axis square
    title('Glauber')
    xlabel('Temperature')
    ylabel('\chi')
    
    
    subplot(2,2,3)
    scatter3(T,PhiMet,PhiSusMet,'MarkerFaceColor',[0 .75 .75])
    xlabel('T')
    ylabel('\Phi')
    zlabel('\chi (\Phi)')
    axis square
    view(90,0)
    
    subplot(2,2,4)
    scatter3(T,PhiGlaub,PhiSusGlaub,'MarkerFaceColor',[0 .75 .75])
    xlabel('T')
    ylabel('\Phi')
    zlabel('\chi (\Phi)')
    axis square
    view(90,0)
    
    suptitle(simFolder{iNet})
    pause()
end
%% Visualize Phi/Sus

% Choose which data to visualize, random or brain networks

data_vis = data_RAND_met;
% data_vis = data_NET_met;

% gifName = 'Ising_Phi_random.gif';

numFiles = length(data_vis);

for i = 1:numFiles % 1:numFiles
    
    J = data_vis(i).J;
    temp = data_vis(i).temp;
    tempPhi = data_vis(i).tempPhi;
    Sus = data_vis(i).Sus;
    Phi = data_vis(i).Phi;
    
    
    figure(1)
    set(gcf, 'Position', [50, 100, 1820, 495]);
    
    subplot(1,3,1)
    imagesc(J)
    title('Structural Connectivity')
    axis square
    
    subplot(1,3,2)
    plot(temp,Sus,'.')
    title('Susceptibility')
    %     axis tight
    
    xlabel(num2str(i))
    
    smoothSus = smooth(temp,Sus,0.1);
    maxSus = max(smoothSus);
    IndC = find(smoothSus == maxSus);
    IndC = IndC(1);
    line([temp(IndC),temp(IndC)],[0,smoothSus(IndC)],'Color','r')
    str = sprintf('$$T_c = %.3g $$',temp(IndC));
    text('Interpreter','latex','Position',[temp(IndC+1) maxSus/10],'String',str,'FontSize',14)
    
    subplot(1,3,3)
    plot(tempPhi,Phi,'.')
    title('Phi')
    %     axis tight
    smoothPhi = smooth(tempPhi,Phi,0.1);
    
    maxPhi = max(smoothPhi);
    IndP = min(find(smoothPhi == maxPhi));
    IndP = IndP(1);
    line([tempPhi(IndP),tempPhi(IndP)],[0,smoothPhi(IndP)],'Color','r')
    str = sprintf('$$T_{max} = %.3g $$',tempPhi(IndP));
    text('Interpreter','latex','Position',[tempPhi(IndP+1) maxPhi/10],'String',str,'FontSize',14)
    
    pause()
    
    %     F(i) = getframe(gcf);
    %
    %     im = frame2im(F(i));
    %     [imind,cm] = rgb2ind(im,256);
    %     if i == 1;
    %         imwrite(imind,cm,gifName,'gif', 'Loopcount',inf);
    %     else
    %         imwrite(imind,cm,gifName,'gif','WriteMode','append');
    %     end
    
end

% fig = figure;
% set(fig, 'Position', [100, 100, 1916, 495]);
% movie(fig,F,2)
%%
clearvars -except data_NET_met data_NET_glaub data_RAND_met data_RAND_glaub Out

normalization = 1; % 1 for normalization, 0 for raw data
markerSize = 15; 
% data = data_RAND_glaub;
data = data_RAND_met;
% data = data_NET_met;
% data = data_NET_glaub;

colors = cbrewer('qual', 'Set2', 8);
fontSize = 20;


temp = data(1).temp;
tempPhi = data(1).tempPhi;

% temp = 1:length(data(1).temp);
% tempPhi = 1:length(data(1).tempPhi);

xt = [0, 1, 2, 3, 4];


Nind{2} = 1:length(data);

i = 2;

catPhi{i} = cat(3,data(Nind{i}).Phi); catPhi{i} = squeeze(catPhi{i});
catSus{i} = cat(3,data(Nind{i}).Sus); catSus{i} = squeeze(catSus{i});
catPhiSus{i} = cat(3,data(Nind{i}).PhiSus); catPhiSus{i} = squeeze(catPhiSus{i});
catM{i} = cat(3,data(Nind{i}).Mag); catM{i} = squeeze(catM{i});

switch normalization
    
    case 0
        
        plotSus = catSus{2};
        stdSus = std(plotSus');
        
        plotPhi = catPhi{2};
        stdPhi = std(plotPhi');
        
        plotPhiSus = catPhiSus{2};
        stdPhiSus = std(plotPhiSus');
        
        plotM = catM{2};
        stdM = std(plotM');
        
    case 1
        
        %%%%%%%% NORMALIZED VERSIONS %%%%%%%%
        
        plotSus = catSus{2};
        plotSus = plotSus./max(plotSus);
        stdSus = std(plotSus');
        
        plotPhi = catPhi{2};
        plotPhi = plotPhi./max(plotPhi);
        stdPhi = std(plotPhi');
        
        plotPhiSus = catPhiSus{2};
        plotPhiSus = plotPhiSus./max(plotPhiSus);
        stdPhiSus = std(plotPhiSus');
        
        plotM = catM{2};
        plotM = plotM./max(plotM);
        stdM = std(plotM');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

figure(1)
set(gcf, 'Position', [300, 50, 1200, 950]);

subplot(2,4,8)
semilogx(tempPhi,stdPhiSus,'.','MarkerSize',markerSize)
title('$\sigma^2_J(\sigma^2_t(\Phi))$','Interpreter','Latex')
box on
grid on
grid minor
axis tight
axis square
set(gca,'FontSize',fontSize)

subplot(2,4,6)
semilogx(tempPhi,stdPhi,'.','MarkerSize',markerSize)
title('$\sigma^2_J(\Phi)$','Interpreter','Latex')
box on
grid on
grid minor
axis tight
axis square
set(gca,'FontSize',fontSize)
set(gca,'GridAlpha', 0.25)

subplot(2,4,4)
semilogx(temp,stdSus,'.','MarkerSize',markerSize)
title('$\sigma^2_J(\chi)$','Interpreter','Latex')
box on
grid on
grid minor
axis tight
axis square
set(gca,'FontSize',fontSize)

subplot(2,4,2)
semilogx(temp,stdM,'.','MarkerSize',markerSize)
title('$\sigma^2_J(M)$','Interpreter','Latex')
box on
grid on
grid minor
axis tight
axis square
set(gca,'FontSize',fontSize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure

subplot(2,4,7)
bl = boundedline(tempPhi, mean(plotPhiSus,2), stdPhiSus,'.-', ...
    'cmap', colors, ...
    'alpha');
title('$\bar{\sigma^2_t}(\Phi)$','Interpreter','Latex')
set(gca,'XScale','log')
box on
grid on
grid minor
axis tight
axis square
set(gca,'FontSize',fontSize)


subplot(2,4,5)
b2 = boundedline(tempPhi, mean(plotPhi,2), stdPhi, '.-', ...
    'cmap', colors, ...
    'alpha');
set(gca,'XScale','log')
title('$\bar{\Phi}$','Interpreter','Latex')
box on
grid on
grid minor
axis tight
axis square
set(gca,'FontSize',fontSize)

subplot(2,4,3)
b3 = boundedline(temp, mean(plotSus,2), stdSus, '.-', ...
    'cmap', colors, ...
    'alpha');
set(gca,'XScale','log')
title('$\bar{\chi}$','Interpreter','Latex')
box on
grid on
grid minor
axis tight
axis square
set(gca,'FontSize',fontSize)

subplot(2,4,1)
b4 = boundedline(temp, mean(plotM,2), stdM, '.-', ...
    'cmap', colors, ...
    'alpha');
set(gca,'XScale','log')
title('$\bar{M}$','Interpreter','Latex')
box on
grid on
grid minor
axis tight
axis square
set(gca,'FontSize',fontSize)

set(bl,'MarkerSize',markerSize)
set(b2,'MarkerSize',markerSize)
set(b3,'MarkerSize',markerSize)
set(b4,'MarkerSize',markerSize)

% if normalization
%     
%     for iSubplot = 1:8
%         subplot(2,4,iSubplot) % normalize yaxis to be between 0 and 1
%         ylim([0 1])
%     end
%     
% end

for iSubplot = 1:8
    
    figure(1)
    subplot(2,4,iSubplot)
    xticks(xt)
    xlim([0 max(xt)])
    
end

figure(1)
l1 = suplabel(['Summary Statistical of ', num2str(length(data)), ' Random Networks'],'t');
% suptitle('Random Networks, Glauber Transitions, n = 28')
l2 = suplabel('Temperature','x',[.1 .12 .84 .84]);

set(l1, 'FontSize',20)
set(l2, 'FontSize',20)
% 
% figure
% plot(tempPhi,plotPhi,'.','MarkerSize',10)
% title('$\Phi$','Interpreter','Latex')
% box on
% grid on
% grid minor
% axis square
% set(gca,'FontSize',fontSize)
% 
% figure
% plot(temp,plotM,'.')
% title('$M$','Interpreter','Latex')
% box on
% grid on
% grid minor
% axis square
% set(gca,'FontSize',fontSize)



% filename = 'figures_paper/rand_net_stats';
% export_fig(filename, '-pdf', '-tiff', '-painters')



%% Averages & Scatter plots (Graph Strenth vs. Phi)
% clearvars -except data data_NET Out
clearvars -except data_NET_met data_NET_glaub data_RAND_met data_RAND_glaub

data = data_RAND_met;
% data = data_RAND_glaub;


temp = data(1).temp;
tempPhi = data(1).tempPhi;
% organize by size N
% Nind{1} = cat(1,data.N) == 5;
% Nind{2} = cat(1,data.N) == 6;
% Nind{3} = cat(1,data.N) == 7;

% organize by batch name
% Nind{1} = strcmp({data.Source} , 'N5');
% Nind{2} = strcmp({data.Source} , 'N5_motif_full_01');
% Nind{2} = strcmp({data.Source} , 'N5_motif_full_newPhi');
% Nind{2} = strcmp({data.Source} , 'N5_motif_full');

Nind{2} = 1:length(data);

colors = cbrewer('qual', 'Set2', 8);

for i = 2
    N = 5; %i + 4;
    
    catPhi{i} = cat(3,data(Nind{i}).Phi); catPhi{i} = squeeze(catPhi{i});
    catSus{i} = cat(3,data(Nind{i}).Sus); catSus{i} = squeeze(catSus{i});
    catPhiSus{i} = cat(3,data(Nind{i}).PhiSus); catPhiSus{i} = squeeze(catPhiSus{i});
    J{i} = cat(3,data(Nind{i}).J);
    %%%% Graph Strength %%%%
    % BEWARE FOR NORMALIZATION!!
    MGS{i} = mean(matrix_squareform(J{i}),1);
    %     normalize = mean(Out.Tmax_RAND,2);
    normalize = 1;
    %%%
    
    MGS{i} = MGS{i}./normalize;
    %%%% # of Edges %%%%
    sparsity{i} = matrix_squareform(J{i});
    sparsity{i}( abs(sparsity{i}) > 0 ) = 1;
    sparsity{i}( ~( abs(sparsity{i}) > 0 ) ) = 0;
    sparsity{i} = sum(sparsity{i},1);
    
    
    
    %%%% Density %%%%
    % sparsity{i} = matrix_squareform(sparsity{i});
    % sparsity{i} = sum(logical(sparsity{i}),1);
    
    %     sparsity{i} = sparsity{i}./nchoosek(N,2);
    
    % Average over simulations that generated significant Phi
    ind = sum(catPhi{i});
    ind = ind > 0.01;
    
    % Plot Sus/Phi average over all random networks that generated non-zero phi
    
    %     figure(i)
    %     plotSus = mean(catSus{i}(:,ind),2)./max(max(mean(catSus{i}(:,ind),2)));
    %     plot(temp,plotSus,'.','MarkerSize',25)
    %     hold on
    %     plotPhi = mean(catPhi{i}(:,ind),2)./max(max(mean(catPhi{i}(:,ind),2)));
    %     plot(temp,plotPhi,'.','MarkerSize',25)
    %     legend({'\chi','\Phi'})
    %     % title(['Average over all N = ',num2str(N),' random networks']);
    %     title(['Average statistics']);
    %     xlabel('Temperature')
    %     set(gca,'FontSize',20)
    
    %%
    %%%%%%%%%%%%%%%%%% AVERAGES %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(i)
%     plotSus = catSus{2}./max(mean(catSus{2},2));    
    plotSus = catSus{2}./max(catSus{2});   
%     stdSus = std(mean(plotSus,2)');
    stdSus = std(plotSus');
    %     plot(temp,mean(plotSus,2),'.','MarkerSize',25)
    %     hold on
%     plotPhi = catPhi{2}./max(mean(catPhi{2},2));    
    plotPhi = catPhi{2}./max(catPhi{2});    
%     stdPhi = std(mean(plotPhi,2)');
    stdPhi = std(plotPhi');
    %     plot(temp,mean(plotPhi,2),'.','MarkerSize',25)
%     plotPhiSus = catPhiSus{2}./max(mean(catPhiSus{2},2));    
    plotPhiSus = catPhiSus{2}./max(catPhiSus{2});  
%     stdPhiSus = std(mean(plotPhiSus,2)');
    stdPhiSus = std(plotPhiSus');
    
%     bl = boundedline(temp, mean(plotSus,2), stdSus, ...
%         tempPhi, mean(plotPhi,2), stdPhi, ...
%         tempPhi, mean(plotPhiSus,2), stdPhiSus, ...
%         'cmap', colors, ...
%         'alpha');
    
        bl = boundedline(tempPhi, mean(plotPhi,2), stdPhi, ...
        'cmap', colors, ...
        'alpha');
    hold on
%     b2 = plot(temp,plotSus,'.','Color', colors(1,:));
%     b3 = plot(tempPhi,plotPhi,'.','Color', colors(2,:));
%     b4 = plot(tempPhi,plotPhiSus,'.','Color',colors(3,:));
    
%     bl = boundedline(temp, plotSus, stdSus,'.', ...
%         temp, plotPhi, stdPhi, '.', ...
%         'cmap', colors, ...
%         'alpha');
    
    set(bl,'LineWidth',10)
    set(b2,'MarkerSize',5)
    set(b3,'MarkerSize',5)
    
    lh = legend(bl);
    legnames = {'\chi : Susceptibility', '\phi : Integrated Information', '\sigma^2(\Phi): std.'};
    for i = 1:length(legnames),
        str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors(i, 1), colors(i, 2), colors(i, 3), legnames{i})];
    end
    lh.String = str;
    lh.FontSize = 15;
    
    %     legend({'\chi','\phi'})
    %     title(['Average over all N = ',num2str(N),' random networks']);
    title(['Random Networks (n = ',num2str(length(data)) ,')']);
    xlabel('Temperature')
    ylabel('au')
    axis tight
%     ylim([0 1.3])
    set(gca,'FontSize',10)
    
    set(gca, ...
        'FontSize'    , 10        , ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.02 .02] , ...
        'XMinorTick'  , 'on'      , ...
        'YMinorTick'  , 'on'      , ...
        'YGrid'       , 'on'      , ...
        'XColor'      , [.3 .3 .3], ...
        'YColor'      , [.3 .3 .3], ...
        'YTick'       , 0:0.25:1.25, ...
        'YTickLabeL'  , [],       ...
        'LineWidth'   , 1         );
    
    % set(gcf, 'PaperPositionMode', 'auto');
    % print(gcf, '-dpdf', 'fig1.pdf');
    % print -depsc2 fig11.eps
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    % %   Sparsity vs. max(Phi) Scatter
    figure(5)
    
%     subplot(1,2,1)
    
    hold on
    plot(MGS{i},max(catPhi{i}),'.','MarkerSize',11)
    title('Mean Graph Strength vs. \Phi_{max}')
    xlabel('Mean Graph Strength')
    ylabel('\Phi_{max}')
    hold off
    set(gca,'FontSize',14)
    axis tight
    axis square
    
%     subplot(1,2,2)
%     
%     hold on
%     plot(sparsity{i},max(catPhi{i}),'.','MarkerSize',11)
%     title('# of Edges vs. \Phi_{max}')
%     xlabel('# of Edges')
%     ylabel('\Phi_{max}')
%     hold off
%     set(gca,'FontSize',14)
%     axis tight
%     axis square
    
    figure(6)
    
%     subplot(1,2,1)
    
    hold on
    plot(MGS{i},mean(catPhi{i}),'.','MarkerSize',11)
    title('Mean Graph Strength vs. \Phi_{mean}')
    xlabel('Mean Graph Strength')
    ylabel('\Phi_{mean}')
    hold off
    set(gca,'FontSize',14)
    axis tight
    axis square
    
%     subplot(1,2,2)
%     
%     hold on
%     plot(sparsity{i},mean(catPhi{i}),'.','MarkerSize',11)
%     title('# of Edges vs. \Phi_{mean}')
%     xlabel('# of Edges')
%     ylabel('\Phi_{max}')
%     hold off
%     set(gca,'FontSize',14)
%     axis tight
%     axis square
    
    % %   Sus vs. Phi Scatter
    
    %     figure(4)
    %     title('Scatter Sus vs. Phi')
    %     ylabel('Phi')
    %     xlabel('Sus')
    %
    %     for ii = length(temp):-1:1
    %         hold all
    %         plot(catSus{1}(ii,ind),catPhi{i}(ii,ind),'.')
    % %         ylim([0 1])
    % %         xlim([0 0.5])
    %         drawnow;
    %     end
    %     hold off
    % %
    %     pause()
end

%% Add Brain Networks to scatter plot (data_net version, legacy)
% Continued from scatter plots above
% RUN paper_figure_maker.m to get data_net structure


% load network data
catPhi{4} = squeeze(cat(3,data_net.Phi));
catSus{4} = squeeze(cat(3,data_net.Sus));
% mean graph strength
MGS{4} = cat(3,data_net.J);
MGS{4} = mean(matrix_squareform(MGS{4}),1);
J = squeeze(cat(3,data_net.J));
% number of edges
sparsity{4} = matrix_squareform(J);
sparsity{4}( abs(sparsity{4}) > 0 ) = 1;
sparsity{4}( ~( abs(sparsity{4}) > 0 ) ) = 0;
sparsity{4} = sum(sparsity{4},1);

% add to random graph scatter plots

figure(5)
% subplot(1,2,1)
hold on
plot(MGS{4},max(catPhi{4}),'o','LineWidth',4)
% 
% subplot(1,2,2)
% hold on
% plot(sparsity{4},max(catPhi{4}),'o','LineWidth',4)

figure(6)
% subplot(1,2,1)
hold on
plot(MGS{4},mean(catPhi{4}),'o','LineWidth',4)

% subplot(1,2,2)
% hold on
% plot(sparsity{4},mean(catPhi{4}),'o','LineWidth',4)

for iNet = 1:9
    figure(5)
%     subplot(1,2,1)
    txt1 = ['\leftarrow ', data_net(iNet).Network ];
    text(MGS{4}(iNet),max(catPhi{4}(:,iNet)),txt1,'FontSize', 14)
    
%     subplot(1,2,2)
%     txt1 = ['\leftarrow ', data_net(iNet).Network ];
%     text(sparsity{4}(iNet),max(catPhi{4}(:,iNet)),txt1,'FontSize', 14)
    
    
    figure(6)
%     subplot(1,2,1)
    txt2 = ['\leftarrow ', data_net(iNet).Network ];
    text(MGS{4}(iNet),mean(catPhi{4}(:,iNet)),txt2,'FontSize',14)
    
%     subplot(1,2,2)
%     txt1 = ['\leftarrow ', data_net(iNet).Network ];
%     text(sparsity{4}(iNet),mean(catPhi{4}(:,iNet)),txt1,'FontSize', 14)
    
end

%
figure(5)
% subplot(1,2,1)
% legend({'N = 5', 'N = 6', 'N = 7', 'Brain Networks'})
legend({'Random motif','Full motif','Brain Networks'})
legend('Location','southeast')
legend('boxoff')
set(gca,'FontSize',14)

% subplot(1,2,2)
% legend({'Random motif','Full motif','Brain Networks'})
% legend('Location','southeast')
% legend('boxoff')
% set(gca,'FontSize',14)

figure(6)
% subplot(1,2,1)
legend({'Random motif','Full motif','Brain Networks'})
legend('Location','southeast')
legend('boxoff')
set(gca,'FontSize',14)

% subplot(1,2,2)
% legend({'Random motif','Full motif','Brain Networks'})
% legend('Location','southeast')
% legend('boxoff')
% set(gca,'FontSize',14)

%% Add Brain Networks to scatter plot (data_NET version)
% Continued from scatter plots above


load('/home/heysoos/Dropbox/Masters Thesis/Codes/Phi/reduced_networks/indices.mat')
simFolder = {'Aud', 'DMN', 'ECN_L', 'ECN_R', 'Salience', 'Sensorimotor', 'VISL', 'VISM', 'VISO'}

clear J
catPhi{4} = [];
catSus{4} = [];
for iNet = 1:length(simFolder)
    Nind = strcmp({data_NET.Network} , simFolder{iNet});
    
    mPhi = mean ( cat(1, data_NET(Nind).Phi), 1);
    mSus = mean ( cat(1, data_NET(Nind).J ), 1);
    
%     figure(99)
%     subplot(1,2,1)
%     plot( cat(1, data_NET(Nind).Phi)', '.')
%     subplot(1,2,2)
%     plot( cat(1, data_NET(Nind).Sus)', '.')
%     pause()
    
    Jt = cat(3,data_NET(Nind).J);
    J(:,:,iNet) = Jt(:,:,1); % no need to avg. just take the first J
    
    catPhi{4}(:,iNet) = mPhi;
    catSus{4}(:,iNet) = mSus;
    
end

% mean graph strength
% BEWARE FOR NORMALIZATION!!
MGS{4} = mean(matrix_squareform(J),1);
% normalize = mean(Out.Tmax_NET,2)';
% set to 1 for default
normalize = 1;
%%%

MGS{4} = MGS{4}./normalize;
% number of edges
sparsity{4} = matrix_squareform(J);
sparsity{4}( abs(sparsity{4}) > 0 ) = 1;
sparsity{4}( ~( abs(sparsity{4}) > 0 ) ) = 0;
sparsity{4} = sum(sparsity{4},1);

% add to random graph scatter plots

figure(5)
% subplot(1,2,1)
hold on
plot(MGS{4},max(catPhi{4}),'o','LineWidth',4)

% subplot(1,2,2)
% hold on
% plot(sparsity{4},max(catPhi{4}),'o','LineWidth',4)

figure(6)
% subplot(1,2,1)
hold on
plot(MGS{4},mean(catPhi{4}),'o','LineWidth',4)

% subplot(1,2,2)
% hold on
% plot(sparsity{4},mean(catPhi{4}),'o','LineWidth',4)

offset = 0;

for iNet = 1:9
    figure(5)
%     subplot(1,2,1)
    txt1 = [simFolder(iNet) ];
    text(MGS{4}(iNet),max(catPhi{4}(:,iNet)) - offset ,txt1,'FontSize', 14)
    
%     subplot(1,2,2)
%     txt1 = ['\leftarrow ', simFolder(iNet) ];
%     text(sparsity{4}(iNet),max(catPhi{4}(:,iNet)),txt1,'FontSize', 14)
    
    
    figure(6)
%     subplot(1,2,1)
%     txt2 = ['\leftarrow ', simFolder(iNet) ];
    txt2 = [simFolder(iNet) ];
    text(MGS{4}(iNet),mean(catPhi{4}(:,iNet)) - offset ,txt2,'FontSize',14)
    
%     subplot(1,2,2)
%     txt1 = ['\leftarrow ', simFolder(iNet) ];
%     text(sparsity{4}(iNet),mean(catPhi{4}(:,iNet)),txt1,'FontSize', 14)
    
end

%
figure(5)
% subplot(1,2,1)
% legend({'N = 5', 'N = 6', 'N = 7', 'Brain Networks'})
legend({'Random Full motif','Brain Networks'})
legend('Location','southeast')
legend('boxoff')
set(gca,'FontSize',14)

% subplot(1,2,2)
% legend({'Random Full motif','Brain Networks'})
% legend('Location','southeast')
% legend('boxoff')
% set(gca,'FontSize',14)

figure(6)
% subplot(1,2,1)
legend({'Random Full motif','Brain Networks'})
legend('Location','southeast')
legend('boxoff')
set(gca,'FontSize',14)

% subplot(1,2,2)
% legend({'Random Full motif','Brain Networks'})
% legend('Location','southeast')
% legend('boxoff')
% set(gca,'FontSize',14)



%% Lines of best fit (graph strength vs. phi_max/mean)
clearvars -except data data_NET Out catPhi MGS

numBatches = 2; % select the data set which you want to make lines of best fit for
% this value is sort of arbitrary, it depends on which cell index for the
% catPhi cell array you are interested in..

% at the moment, 2 is for the random data set, 4 is for the networks
for i = 2
    holderPhiMax{i} = max(catPhi{i});
    holderPhiMean{i} = mean(catPhi{i});
end

randPhiMax = cat(2,holderPhiMax{numBatches});
randPhiMean = cat(2,holderPhiMean{numBatches});

randSpars = cat(2,MGS{numBatches});

fitMax = polyfit(randSpars,randPhiMax,1);
fitMean = polyfit(randSpars,randPhiMean,1);

yMax = polyval(fitMax,randSpars);
yMean = polyval(fitMean,randSpars);

% compute residuals

resMax = randPhiMax - yMax;
SSresMax = sum(resMax.^2);
SStotalMax = (length(randPhiMax) - 1) * var(randPhiMax);

% R^2_Max
rsqMax = 1 - SSresMax/SStotalMax;

resMean = randPhiMean - yMean;
SSresMean = sum(resMean.^2);
SStotalMean = (length(randPhiMean) - 1) * var(randPhiMean);

% R^2()_Mean
rsqMean = 1 - SSresMean/SStotalMean;

% graph lines of best fit

x = 0:0.01:0.9;

figure(5)
% subplot(1,2,1)
plot(x,polyval(fitMax,x),'--','LineWidth',3,'Color','black')
txt1 = ['R^2 = ', num2str(rsqMax)];
text(0.5,0.95,txt1,'FontSize',14,'Units','normalized','HorizontalAlignment','center')

figure(6)
% subplot(1,2,1)
plot(x,polyval(fitMean,x),'--','LineWidth',3,'Color','black')
txt2 = ['R^2 = ', num2str(rsqMean)];
text(0.5,0.95,txt2,'FontSize',14,'Units','normalized','HorizontalAlignment','center')

%% Correlation distance of Random Sims/Brain Networks to Empirical Correlations
% To choose between Random Sims and Brain Networks, comment the appropriate
% comment blocks labelled.

clearvars -except data_NET_met data_NET_glaub data_RAND_met data_RAND_glaub Out

homePhi = '../'; % Should lead to main 'Phi' directory
Corr_FMRI = load([homePhi, '/connectome_data/mean_struct_corr.mat']);
Corr_FMRI = Corr_FMRI.mean_corr_fc;
load([homePhi, 'reduced_networks/indices_84.mat'])

Networks = {'Aud', 'DMN', 'ECN_L', 'ECN_R', 'Salience', 'Sensorimotor', 'VISL', 'VISM', 'VISO'};

smoothSpan = 0.1; % the smoothing span with which the maximas/minimas will be found with

% CORRELATION DISTANCES

% COMMENT (Brain Network Simulations) [
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Brain Network Simulations

data_NET = data_NET_met;

Corr_DTI = cat(4,data_NET(:).Corr_DTI);
Sus = cat(1,data_NET(:).Sus);
Phi = cat(1,data_NET(:).Phi);
PhiSus = cat(1,data_NET(:).PhiSus);

temp = cat(1,data_NET(:).temp);
temp = temp(1,:);

tempPhi = cat(1,data_NET(:).tempPhi);
tempPhi = tempPhi(1,:);

% calculate the distances between all random sims with the 9 empirical
% correlations from the networks

ks2stat_NET = [];
PD_NET = [];
Tmin_NET = [];
Nind = [];
pKS_NET = [];
for iNet = 1:length(Networks)
    iNet
    Nind = strfind({data_NET.Network} , Networks{iNet});
    Nind = find(~cellfun(@isempty,Nind));
    
    Ind = indices.(Networks{iNet});
    iterCounter = 0;
    for iter = 1:length(Nind) % sims loop ( SKIPPING BY 3s!!! )
        
        % smooth data to find min/max consistantly
        smoothSus = smooth(temp,Sus(Nind(iter),:),smoothSpan);
        [~, IndD] = max(smoothSus);
        
        smoothPhi = smooth(tempPhi,Phi(Nind(iter),:),smoothSpan);
        [~, IndP] = max(smoothPhi);
        
        smoothPhiSus = smooth(tempPhi,PhiSus(iter,:),smoothSpan);
        [~, IndPS] = max(smoothPhiSus);
        %         [~, IndPS] = max(PhiSus(Nind(iter),:));
        
        %%%%%%%%%%%% DEBUG ONLY %%%%%%%%%%
        
        %         subplot(1,2,1)
        %         plot(temp,Sus(Nind(iter),:),'.',temp,smoothSus)
        %         title('\chi')
        %
        %         subplot(1,2,2)
        %         plot(tempPhi,Phi(Nind(iter),:),'.',tempPhi,smoothPhi)
        %         title('\Phi')
        %
        %         suptitle(Networks(iNet))
        %         pause()
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Tc_NET(iNet,iter) = temp(IndD);
        Tmax_NET(iNet,iter) = temp(IndP);
        Tcphi_NET(iNet,iter) = temp(IndPS);
        
        iterCounter = iterCounter + 1;
        [ks2stat_NET(iNet,iterCounter,:), pKS_NET(iNet,iterCounter,:), PD_NET(iNet,iterCounter,:)] = ks_test_dist(Corr_DTI(:,:,:,Nind(iter)),Corr_FMRI(Ind,Ind));
        smoothKS = smooth(temp,ks2stat_NET(iNet,iterCounter,:),0.1);
        smoothPD = smooth(temp,PD_NET(iNet,iterCounter,:),0.1);
        
        [~,InDM] = min(PD_NET(iNet,iter,:));
        Tmin_NET(iNet,iter) = temp(InDM);
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMMENT (Brain Network Simulations) ]

% COMMENT (Random Network Simulations ) [
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Random Network Simulations

data = data_RAND_met;

Nind = [];
Nind{1} = cat(1,data.N) == 5;

Corr_DTI = cat(4,data(Nind{1}).Corr_DTI);
Sus = cat(1,data(Nind{1}).Sus);
Phi = cat(1,data(Nind{1}).Phi);
PhiSus = cat(1,data(Nind{1}).PhiSus);

temp = cat(1,data(:).temp);
temp = temp(1,:);

tempPhi = cat(1,data(:).tempPhi);
tempPhi = tempPhi(1,:);

N = size(Corr_DTI,1);
nSims = size(Corr_DTI,4);
nFiles = nSims*size(perms(1:N),1); 
nNetworks = length(Networks);
nDataLen = size(Corr_DTI,3);

ks2stat_RAND = zeros(nNetworks, nFiles, nDataLen);
PD_RAND = zeros(nNetworks, nFiles, nDataLen);
Tmin_RAND = zeros(nNetworks, nFiles, nDataLen);
pKS_RAND = zeros(nNetworks, nFiles, nDataLen);

JPermsInd = perms(1:5); % permuations to cycle through for each random sim
TminCounter = 0;

for iter = 1:nSims 
    % smooth data to find min/max consistantly
    smoothSus = smooth(temp,Sus(iter,:),smoothSpan);
    [~, IndD] = max(smoothSus);
    
    smoothPhi = smooth(tempPhi,Phi(iter,:),smoothSpan);
    [~, IndP] = max(smoothPhi);
    
    smoothPhiSus = smooth(tempPhi,PhiSus(iter,:),smoothSpan);
    [~, IndPS] = max(smoothPhiSus);
    %     [~, IndPS] = max(PhiSus(iter,:));
    
    
    %%%%%%%%%%%% DEBUG ONLY %%%%%%%%%%
    
    %     subplot(1,2,1)
    %     plot(temp,Sus(iter,:),temp,smoothSus)
    %     title('\chi')
    %
    %     subplot(1,2,2)
    %     plot(tempPhi,Phi,tempPhi,smoothPhi)
    %     title('\Phi')
    %
    %     pause()
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Tc_RAND(iter) = temp(IndD);
    Tmax_RAND(iter) = temp(IndP);
    Tcphi_RAND(iter) = temp(IndPS);
    
    for iPerms = 1:size(JPermsInd,1)
        TminCounter = TminCounter + 1;
        
        txtStr = [num2str(iter), ' - ', num2str(iPerms)];
        display(txtStr)
        
        Corr_DTI_temp = squeeze(squeeze(Corr_DTI(:,:,:,iter)));
        Corr_DTI_temp = Corr_DTI_temp(JPermsInd(iPerms,:),JPermsInd(iPerms,:),:);
        
        for iNet = 1:length(Networks)
            
            Ind = indices.(Networks{iNet});
            
            [ks2stat_RAND(iNet,TminCounter,:), pKS_RAND(iNet,TminCounter,:), PD_RAND(iNet,TminCounter,:)] = ks_test_dist(Corr_DTI_temp,Corr_FMRI(Ind,Ind));
            
            [~,InDM] = min(PD_RAND(iNet,TminCounter,:));
            Tmin_RAND(iNet,TminCounter) = temp(InDM);
            %             smoothks2stat = smooth(temp,ks2stat,0.1); % needed to locate
            %             % minimums smoothly
            
            
            
            
        end
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMMENT (Random Network Simulations ) ]

% save output in data struct
Out.PD_NET = PD_NET;
Out.PD_RAND = PD_RAND;

Out.ks2stat_NET = ks2stat_NET;
Out.ks2stat_RAND = ks2stat_RAND;

Out.Tmin_NET = Tmin_NET;
Out.Tmin_RAND = Tmin_RAND;

Out.Tc_NET = Tc_NET;
Out.Tc_RAND = Tc_RAND;

Out.Tmax_NET = Tmax_NET;
Out.Tmax_RAND = Tmax_RAND;

Out.Tcphi_NET = Tcphi_NET;
Out.Tcphi_RAND = Tcphi_RAND;


clearvars -except data_NET_met data_NET_glaub data_RAND_met data_RAND_glaub Out

%% Visualize Correlation Distances to respective networks(choose data)
%%---------- SET DATA TO BE VISUALIZED (IMPORTANT!) ----------%%
% set one consistant temperature reference frame. IMPORTANT DETAIL!
% BE CAREFUL!

homePhi = '../'; % Should lead to main 'Phi' directory
Corr_FMRI = load([homePhi, '/connectome_data/mean_struct_corr.mat']);
Corr_FMRI = Corr_FMRI.mean_corr_fc;
load([homePhi, 'reduced_networks/indices_84.mat'])

Networks = {'Aud', 'DMN', 'ECN_L', 'ECN_R', 'Salience', 'Sensorimotor', 'VISL', 'VISM', 'VISO'};

temp = cat(1,data_NET_met(:).temp);
temp = temp(1,:);


ks2stat = Out.ks2stat_NET;
PD = Out.PD_NET;

% ks2stat = Out.ks2stat_NET;
% PD = Out.PD_NET;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% KSTEST2
% plot avged distance vs. temp for all networks
% average across random sims
meanKS = squeeze(mean(ks2stat,2));

for i = 1:9
    figure(1)
    subplot(3,3,i)
    set(gca,'FontSize',14)
    plot(temp,meanKS(i,:),'.')
    title(Networks{i})
    axis tight
    %     ylim([0 1])
    
    [minDist, IndD] = min(meanKS(i,:));
    line([temp(IndD),temp(IndD)],[0,meanKS(i,IndD)],'Color','r')
    str = sprintf('$$T_{min} = %.3g $$',temp(IndD));
    text('Interpreter','latex','Position',[temp(IndD+1) minDist/2],'String',str,'FontSize',15)
end
suptitle('Correlation ks-test vs. Temp')

% PDIST
% plot avged distance vs. temp for all networks
meanPD = squeeze(mean(PD,2));
for i = 1:9
    figure(2)
    subplot(3,3,i)
    set(gca,'FontSize',14)
    plot(temp,meanPD(i,:),'.')
    title(Networks{i})
    axis tight
    %     ylim([0 1])
    
    [minDist, IndD] = min(meanPD(i,:));
    line([temp(IndD),temp(IndD)],[0,meanPD(i,IndD)],'Color','r')
    str = sprintf('$$T_{min} = %.3g $$',temp(IndD));
    text('Interpreter','latex','Position',[temp(IndD+1) minDist/2],'String',str,'FontSize',15)
end
suptitle('Correlation Pdist vs. Temp')

% COMMENT (Visualize one simulation at a time) [
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ks2stat = ks2stat_NET;
% temp = temp(1,:);
% for iter = 1:size(ks2stat,2);
%     KS = squeeze(ks2stat(:,iter,:));
%     for i = 1:9
%         figure(7)
%         subplot(3,3,i)
%         plot(temp,KS(i,:),'.')
%         title(Networks{i})
%         axis tight
%         ylim([0 1])
%
%         [minDist, IndD] = min(KS(i,:));
%         line([temp(IndD),temp(IndD)],[0,KS(i,IndD)],'Color','r')
%         str = sprintf('$$T_{min} = %.3g $$',temp(IndD));
%         text('Interpreter','latex','Position',[temp(IndD+1) minDist/2],'String',str,'FontSize',15)
%     end
%     suptitle(['#: ',num2str(iter) ])
%     pause()
%
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMMENT (Visualize one simulation at a time) ]

clearvars -except data_NET_met data_NET_glaub data_RAND_met data_RAND_glaub Out

%% Histograms of distances (run distance calculations first, above, to get variable PD)

fontSize = 16;

temp = cat(1,data_NET_met(:).temp);
temp = temp(1,:);

PD_RAND = Out.PD_RAND;
Tc_RAND = Out.Tc_RAND;

PD_NET = Out.PD_NET;
Tc_NET = Out.Tc_NET;

Networks = {'Aud', 'DMN', 'ECN_L', 'ECN_R', 'Salience', 'Sensorimotor', 'VISL', 'VISM', 'VISO'};

binWidth = 0.02;
binWidth_NET = 0.02;

distMin_RAND = [];
distMin_NET = [];

numSims = length(Out.Tc_RAND);
for iNet = 1:9
    % random sims
    for iSim = 1:size(PD_RAND,2)
        
        xRand = squeeze(PD_RAND(iNet,iSim,:));
        
        % minimum distances
        [distMin_RAND(iNet,iSim), I_RAND] = min(xRand);
        
        % critical distances
        simIndex = mod(iSim,numSims);
        if simIndex == 0
            simIndex = numSims;
        end
        TcInd = find(temp == Tc_RAND( simIndex ));
        distc_RAND(iNet,iSim) = xRand(TcInd);
        
        %         plot(temp,x)
        %         line([temp(I),temp(I)],[0,distMin(iNet,iSim)],'Color','r')
        %         pause()
        
    end
    % network sims
    for iSim = 1:size(PD_NET,2)
        iSim
        xNET = squeeze(PD_NET(iNet,iSim,:));
        
        % minimum distances
        [distMin_NET(iNet,iSim), I_NET] = min(xNET);
        
%         pause()
        
        %         % critical distances
        %         TcInd = find(temp == Tc_NET(iNet,iSim));
        %         distc_NET(iNet,iSim) = xNET(TcInd);
        %
        %         % max phi distances
        %         TmaxInd = find(temp == Tmax_NET(iNet,iSim));
        %         distMax_NET(iNet,iSim) = xNET(TmaxInd);
    end
    
    % minimum distance distributions
    figure(1)
    subplot(3,3,iNet)
    set(gca,'FontSize',fontSize)
    hold all
    dRAND = distMin_RAND(iNet,:);
    h1 = histogram(dRAND);
    h1.BinWidth = binWidth;
    h1.Normalization = 'probability';
    
    yLimRand = max(h1.BinCounts./sum(h1.BinCounts));
    
    dNET = distMin_NET(iNet,:);
    h2 = histogram(dNET);
    h2.BinWidth = binWidth_NET;
    h2.Normalization = 'probability';
    grid on
    grid minor
    
    p = [];
    for i = 1:length(dNET)
        p(i) = sum(dNET(i)>dRAND)/length(dRAND);
    end
    p = mean(p);
    
    pstr = sprintf(' p = %1.3f', p);
    
    title([Networks{iNet}, pstr])
    
    if iNet == 3
        legend({'Random full networks','Brain networks'})
    end
    
    xlim([0 2])
    ylim([0 yLimRand])
    
    %     % critical distance distributions
    %
    %     figure(2)
    %     subplot(3,3,iNet)
    %     set(gca,'FontSize',14)
    %     hold all
    %     h1 = histogram(distc_RAND(iNet,:));
    %     h1.BinWidth = binWidth;
    %     h1.Normalization = 'probability';
    %
    %     h2 = histogram(distc_NET(iNet,:));
    %     h2.BinWidth = binWidth_NET;
    %     h2.Normalization = 'probability';
    %
    %     title(Networks{iNet})
    %     legend({'Random full networks','Brain networks'})
    %
    %     xlim([0 2])
    %
    %     % maximum distance distributions
    %
    %     figure(3)
    %     subplot(3,3,iNet)
    %     set(gca,'FontSize',14)
    %     hold all
    %     h1 = histogram(distMax_RAND(iNet,:));
    %     h1.BinWidth = binWidth;
    %     h1.Normalization = 'probability';
    %
    %     h2 = histogram(distMax_NET(iNet,:));
    %     h2.BinWidth = binWidth_NET;
    %     h2.Normalization = 'probability';
    %
    %     title(Networks{iNet})
    %     legend({'Random full networks','Brain networks'})
    %
    %     xlim([0 2])
    
    % differences in distances
    
    %     % dist(Tmax) - dist(Tmin)
    %     figure(4)
    %     subplot(3,3,iNet)
    %     set(gca,'FontSize',14)
    %     hold all
    %     h1 = histogram(distMax_RAND(iNet,:) - distMin_RAND(iNet,:));
    %     h1.BinWidth = binWidth;
    %     h1.Normalization = 'probability';
    %
    %     h2 = histogram(distMax_NET(iNet,:) - distMin_NET(iNet,:));
    %     h2.BinWidth = binWidth_NET;
    %     h2.Normalization = 'probability';
    %
    %     title(Networks{iNet})
    %     legend({'Random full networks','Brain networks'})
    %
    %     % dist(Tmax) - dist(Tc)
    %     figure(5)
    %     subplot(3,3,iNet)
    %     set(gca,'FontSize',14)
    %     hold all
    %     h1 = histogram(distMax_RAND(iNet,:) - distc_RAND(iNet,:));
    %     h1.BinWidth = binWidth;
    %     h1.Normalization = 'probability';
    %
    %     h2 = histogram(distMax_NET(iNet,:) - distc_NET(iNet,:));
    %     h2.BinWidth = binWidth_NET;
    %     h2.Normalization = 'probability';
    %
    %     title(Networks{iNet})
    %     legend({'Random full networks','Brain networks'})
    %
    %     % dist(Tc) - dist(Tmin)
    %     figure(6)
    %     subplot(3,3,iNet)
    %     set(gca,'FontSize',14)
    %     hold all
    %     h1 = histogram(distc_RAND(iNet,:) - distMin_RAND(iNet,:));
    %     h1.BinWidth = binWidth;
    %     h1.Normalization = 'probability';
    %
    %     h2 = histogram(distc_NET(iNet,:) - distMin_NET(iNet,:));
    %     h2.BinWidth = binWidth_NET;
    %     h2.Normalization = 'probability';
    %
    %     title(Networks{iNet})
    %     legend({'Random full networks','Brain networks'})
    
    
    % temperature differences
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %     % Tmax - Tmin
    %     figure(7)
    %     subplot(3,3,iNet)
    %     set(gca,'FontSize',14)
    %     hold all
    %
    %     a = Out.Tmax_RAND - Out.Tmin_RAND(iNet,:);
    %     h1 = histogram(a);
    %     h1.BinWidth = binWidth;
    %     h1.Normalization = 'probability';
    %
    %     a = Out.Tmax_NET(iNet,:) - Out.Tmin_NET(iNet,:);
    %     h2 = histogram(a);
    %     h2.BinWidth = binWidth_NET;
    %     h2.Normalization = 'probability';
    %
    %     title(Networks{iNet})
    %     legend({'Random full networks','Brain networks'})
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %     % Tmax - Tc
    %         figure(8)
    %         subplot(3,3,iNet)
    %         set(gca,'FontSize',14)
    %         hold all
    %
    %         a = Out.Tmax_RAND - Out.Tc_RAND;
    %         h1 = histogram(a);
    %         h1.BinWidth = binWidth;
    %         h1.Normalization = 'probability';
    %
    %         a = Out.Tmax_NET(iNet,:) - Out.Tc_NET(iNet,:);
    %         h2 = histogram(a);
    %         h2.BinWidth = binWidth_NET;
    %         h2.Normalization = 'probability';
    %
    %         title(Networks{iNet})
    %         legend({'Random full networks','Brain networks'})
    
    %     Tc - Tmin
%     figure(9)
%     subplot(3,3,iNet)
%     set(gca,'FontSize',14)
%     hold all
%     
%     a = Out.Tc_RAND - Out.Tmin_RAND(iNet,:);
%     h1 = histogram(a);
%     h1.BinWidth = binWidth;
%     h1.Normalization = 'probability';
%     
%     a = Out.Tc_NET(iNet,:) - Out.Tmin_NET(iNet,:);
%     h2 = histogram(a);
%     h2.BinWidth = binWidth_NET;
%     h2.Normalization = 'probability';
%     
%     title(Networks{iNet})
%     legend({'Random full networks','Brain networks'})
%     xlim([-1.5 1.5])
%     
%     %     Tc - Tc_phi
%     figure(11)
%     subplot(3,3,iNet)
%     set(gca,'FontSize',14)
%     hold all
%     
%     a = Out.Tc_RAND - Out.Tmin_RAND(iNet,:);
%     h1 = histogram(a);
%     h1.BinWidth = binWidth;
%     h1.Normalization = 'probability';
%     
%     a = Out.Tc_NET(iNet,:) - Out.Tmin_NET(iNet,:);
%     h2 = histogram(a);
%     h2.BinWidth = binWidth_NET;
%     h2.Normalization = 'probability';
%     
%     title(Networks{iNet})
%     legend({'Random full networks','Brain networks'})
%     xlim([-1.5 1.5])
end
figure(1)
l1 = suplabel('Minimum distance','x',[.1 .10 .84 .84]);
l2 = suplabel('Probability','y',[.12 .1 .84 .84]);
l3 = suplabel('Minimum distance distributions','t',[.1 .12 .84 .84]);

set(l1, 'FontSize',20)
set(l2, 'FontSize',20)
set(l3, 'FontSize',20)

% figure(2)
% suptitle('Critical distance distributions')
%
% figure(3)
% suptitle('\Phi_{max} distance distributions')
%
% figure(4)
% suptitle('d(T_{max}) - d(T_{min})')
%
% figure(5)
% suptitle('d(T_{max}) - d(T_{c})')
%
% figure(6)
% suptitle('d(T_{c}) - d(T_{min})')

% figure(7)
% suptitle('T_{max} - T_{min}')

% figure(8)
% suptitle('T_{max} - T_{c}')
%
% figure(9)
% l1 = suplabel('Difference','x',[.1 .10 .84 .84]);
% l2 = suplabel('Probability','y',[.12 .1 .84 .84]);
% l3 = suplabel('T_{c} - T_{min}','t',[.1 .12 .84 .84]);
% 
% set(l1, 'FontSize',16)
% set(l2, 'FontSize',16)
% set(l3, 'FontSize',16)

% filename = 'figures_paper/min_dist_hist';
% export_fig(filename, '-pdf', '-tiff', '-painters')

clearvars -except data_NET_met data_NET_glaub data_RAND_met data_RAND_glaub Out

%% Empirical/Simulated Correlations & Connectivity (per iteration)

homePhi = '../'; % Should lead to main 'Phi' directory
Corr_FMRI = load([homePhi, '/connectome_data/mean_struct_corr.mat']);
Corr_FMRI = Corr_FMRI.mean_corr_fc;
load([homePhi, 'reduced_networks/indices_84.mat'])

Networks = {'Aud', 'DMN', 'ECN_L', 'ECN_R', 'Salience', 'Sensorimotor', 'VISL', 'VISM', 'VISO'};

numS = 1;

% subjectNum = 1:length(data_NET)./9;

for iNet = 1:length(Networks)
    
    Nind = strcmp({data_NET.Network} , Networks(iNet));
    data_1NET = data_NET(Nind);
    
    temp = data_1NET(1).temp; % just grab the first since all identical
    
    J = data_1NET(1).J; % just grab the first J since they're all identical
    
    Corr_DTI = cat(4,data_1NET(numS).Corr_DTI); % grab 1 simulations's Corr_DTI
    Tmin = Out.Tmin_NET;
    Tmin = Tmin(iNet,numS);
    TminI = find(temp==Tmin);
    Corr_DTI = Corr_DTI(:,:,TminI);
    Corr_DTI = squeeze(Corr_DTI);
    
    I = indices.(Networks{iNet});
    Corr_sFMRI = Corr_FMRI(I,I); % generate the relevant network Corr_FMRI
    
    
    figure(1)
    
    subplot(1,3,1)
    imagesc(J)
    axis square
    title('Conncectivity')
    
    subplot(1,3,2)
    imagesc(Corr_DTI)
    axis square
    title('Simulated Correlations (Tmin)')
    xlabel(Tmin)
    
    subplot(1,3,3)
    imagesc(Corr_sFMRI)
    axis square
    title('Empirical Correlations')
    
    suptitle(Networks{iNet})
    pause()
end

clearvars -except data_NET_met data_NET_glaub data_RAND_met data_RAND_glaub Out

%% Empirical/Simulated Correlations & Connectivity (averaged over iterations)

homePhi = '../'; % Should lead to main 'Phi' directory

Corr_FMRI = load([homePhi, '/connectome_data/mean_struct_corr.mat']);
Corr_FMRI = Corr_FMRI.mean_corr_fc;

load([homePhi, 'reduced_networks/indices_84.mat'])

Networks = {'Aud', 'DMN', 'ECN_L', 'ECN_R', 'Salience', 'Sensorimotor', 'VISL', 'VISM', 'VISO'};

data_NET = data_NET_met;
% numS = 1:length(data_NET);
numS = 1;

for iNet = 1:length(Networks)
    
    Nind = strcmp({data_NET.Network} , Networks(iNet));
    data_1NET = data_NET(Nind);
    
    temp = data_1NET(1).temp; % just grab the first since all identical
    
    J = data_1NET(1).J; % just grab the first J since they're all identical
    J(1:length(J)+1:end) = NaN;
    
    Corr_DTI = cat(4,data_1NET(numS).Corr_DTI); % grab 1 simulations's Corr_DTI
    Tmin = Out.Tmin_NET;
    Tmin = Tmin(iNet,:);
    for iS = 1:length(Tmin)
        TminI(iS) = find(temp==Tmin(iS));
        Corr_pDTI(:,:,:,iS) = Corr_DTI(:,:,TminI(iS),iS);
    end
    Corr_pDTI = squeeze(mean(Corr_pDTI,4));
    Corr_pDTI(1:length(Corr_pDTI)+1:end) = NaN;
    
    I = indices.(Networks{iNet});
    Corr_sFMRI = Corr_FMRI(I,I); % generate the relevant network Corr_FMRI
    Corr_sFMRI(1:length(Corr_sFMRI)+1:end) = NaN;
    
    cmin = min([ min(min(J)) min(min(Corr_pDTI)) min(min(Corr_sFMRI)) ]);
    cmax = max([ max(max(J)) max(max(Corr_pDTI)) max(max(Corr_sFMRI)) ]);
    
    figure('units','normalized','position',[.1 .1 .7 .4])
    
    subplot(1,3,1)
    imagesc(J)
    axis square
    title('Conncectivity (DTI)')
    box off
    xticks([])
    yticks([])
    grid on
    caxis([cmin cmax])
    colormap(linspecer)
%     colorbar
    
    subplot(1,3,2)
    imagesc(Corr_pDTI)
    axis square
    title('Mean Ising Correlations (Tmin)')
    txt1 = sprintf('mean T_{min} = %1.2f', mean(Tmin));
    xlabel(txt1)
    box off
    xticks([])
    yticks([])
    grid on
    caxis([cmin cmax])
    colormap(linspecer)
%     colorbar
    
    subplot(1,3,3)
    imagesc(Corr_sFMRI)
    axis square
    title('Empirical Correlations (fMRI)')
    box off
    xticks([])
    yticks([])
    grid on
    caxis([cmin cmax])
    colormap(linspecer)
%     colorbar
    
    suptitle(Networks{iNet})
    pause()
end

clearvars -except data_NET_met data_NET_glaub data_RAND_met data_RAND_glaub Out

%% Temperature scatter plots
Networks = {'Aud', 'DMN', 'ECN_L', 'ECN_R', 'Salience', 'Sensorimotor', 'VISL', 'VISM', 'VISO'};

%%%%%%%%%%%%%%%%%%%%% Tc vs Tmax %%%%%%%%%%%%%%%%%%%%%
figure()
subplot(1,3,1)

Tcr = Out.Tc_RAND;
Tmr = Out.Tmax_RAND;

plot(Tcr,Tmr,'.','MarkerSize',11)
title('T_c vs. T_{max}')
xlabel('T_c')
ylabel('T_{max}')
hold off
set(gca,'FontSize',14)
axis tight
axis square

Tcn = mean(Out.Tc_NET,2);
Tmn = mean(Out.Tmax_NET,2);

hold on

plot(Tcn,Tmn,'.','MarkerSize',11)

fit = polyfit(Tcr,Tmr,1);

y = polyval(fit,Tcr);

% compute residuals

res = Tmr - y;
SSres = sum(res.^2);
SStotal = (length(Tmr) - 1) * var(Tmr);

% R^2_Max
rsq = 1 - SSres/SStotal;

% graph lines of best fit

x = 0:0.01:3;

% plot(x,polyval(fit,x),'--','LineWidth',3,'Color','black')
% txt1 = ['R^2 = ', num2str(rsq)];
% text(0.5,0.95,txt1,'FontSize',14,'Units','normalized','HorizontalAlignment','center')

rho = corr(Tcr',Tmr','Type','Spearman');

% txt2 = ['\rho = ', num2str(rho)];
% text(0.5,0.85,txt2,'FontSize',14,'Units','normalized','HorizontalAlignment','center')

txt3 = sprintf('R^2 = %0.2f\n \\rho = %0.2f', rsq, rho);
text(0.5,0.95,txt3,'FontSize',14,'Units','normalized','HorizontalAlignment','center')

%%%%%%%%%%%%%%%%%%%%% Tc vs. Tmin %%%%%%%%%%%%%%%%%%%%%
subplot(1,3,2)

Tcr = Out.Tc_RAND;
Tmr = Out.Tmin_RAND;

plot(Tcr,Tmr,'.','MarkerSize',11)
title('T_c vs. T_{min}')
xlabel('T_c')
ylabel('T_{min}')
hold off
set(gca,'FontSize',14)
axis tight
axis square

Tcn = mean(Out.Tc_NET,2);
Tmn = mean(Out.Tmin_NET,2);

hold on

plot(Tcn,Tmn,'.','MarkerSize',11)
for iNet = 1:9
    
    fit = polyfit(Tcr,Tmr(iNet,:),1);
    
    y = polyval(fit,Tcr);
    
    % compute residuals
    
    res = Tmr - y;
    SSres = sum(res.^2);
    SStotal = (length(Tmr) - 1) * var(Tmr);
    
    % R^2_Max
    rsq = 1 - SSres/SStotal;
    
    % graph lines of best fit
    
    x = 0:0.01:3;
    
    plot(x,polyval(fit,x),'--','LineWidth',1,'Color',[0.5 0.5 0.5])
    % txt1 = ['R^2 = ', num2str(rsq)];
    % text(0.5,0.95,txt1,'FontSize',14,'Units','normalized','HorizontalAlignment','center')
    
end

legend('Location','northwest')
legend('boxoff')
legend(Networks)
set(gca,'FontSize',14)

%%%%%%%%%%%%%%%%%%%%% Tmax vs. Tmin %%%%%%%%%%%%%%%%%%%%%
subplot(1,3,3)

Tmar = Out.Tmax_RAND;
Tmr = Out.Tmin_RAND;

plot(Tmar,Tmr,'.','MarkerSize',11)
title('T_{max} vs. T_{min}')
xlabel('T_{max}')
ylabel('T_{min}')
hold off
set(gca,'FontSize',14)
axis tight
axis square

Tman = mean(Out.Tmax_NET,2);
Tmn = mean(Out.Tmin_NET,2);

hold on

plot(Tman,Tmn,'.','MarkerSize',11)
for iNet = 1:9
    
    fit = polyfit(Tmar,Tmr(iNet,:),1);
    
    y = polyval(fit,Tmar);
    
    % compute residuals
    
    res = Tmr - y;
    SSres = sum(res.^2);
    SStotal = (length(Tmr) - 1) * var(Tmr);
    
    % R^2_Max
    rsq = 1 - SSres/SStotal;
    
    % graph lines of best fit
    
    x = 0:0.01:3;
    
    plot(x,polyval(fit,x),'--','LineWidth',1,'Color',[0.5 0.5 0.5])
    % txt1 = ['R^2 = ', num2str(rsq)];
    % text(0.5,0.95,txt1,'FontSize',14,'Units','normalized','HorizontalAlignment','center')
    
end

legend('Location','northwest')
legend('boxoff')
legend(Networks)
set(gca,'FontSize',14)


suptitle('Temperature scatter plots')

clearvars -except data_NET_met data_NET_glaub data_RAND_met data_RAND_glaub Out

%% ttest temperatures
Networks = {'Aud', 'DMN', 'ECN_L', 'ECN_R', 'Salience', 'Sensorimotor', 'VISL', 'VISM', 'VISO'};

binWidth = 0.1;

Tc = Out.Tc_RAND;
Tmax = Out.Tmax_RAND;

% IndF = Tc < 0.1;
% IndF2 = Tmax < 0.1;

[h,p] = ttest(Tc,Tmax)
[h,p] = ttest(Tc./Tmax-1)

figure()
subplot(1,2,1)
histogram(Tc)
hold on
histogram(Tmax)
title('Compare Tc and Tmax Distributions')
legend({'Tc','Tmax'})
axis square

subplot(1,2,2)
histogram(Tc./Tmax,200)
title('Tc/Tmax ratio')
% xlim([0 10])
axis square

[h,p] = kstest2(Tc,Tmax)

figure()
for iNet = 1:length(Networks)
    
    subplot(3,3,iNet)
    
    Tmin = Out.Tmin_RAND(iNet,:);
    h1 = histogram(Tc);
    hold all
    h2 = histogram(Tmax);
    h3 = histogram(Tmin);
    h1.BinWidth = binWidth;
    h2.BinWidth = binWidth;
    h3.BinWidth = binWidth;
    legend({'Tc','Tmax','Tmin'})
    
    title(Networks{iNet})
    
end

suptitle('Compare T_c & T_{max} with T_{min}s')

clearvars -except data_NET_met data_NET_glaub data_RAND_met data_RAND_glaub Out

%% Avg. Sus, Phi vs. T Figures
homePhi = '../'; % Should lead to main 'Phi' directory
load([homePhi, 'reduced_networks/indices_84.mat'])

Networks = {'Aud', 'DMN', 'ECN_L', 'ECN_R', 'Salience', 'Sensorimotor', 'VISL', 'VISM', 'VISO'};

data_NET = data_NET_met;

numS = 1:length(data_NET)/9;
C = colormap;
colors = cbrewer('qual','Set2',8);
mkSize = 10;
fontSize = 16;
lWidth = 3;
start = 1;

% plotCase = 'ener';
plotCase = 'phi';
% plotCase = 'sus/spec';
xt = 0:4;


drawLine = true;

drawPhiSus = true;
drawPhiSusLine = false;

drawTcLine = false;

for iNet = 1:length(Networks)
    
    Nind = strcmp({data_NET.Network} , Networks(iNet));
    data_1NET = data_NET(Nind);
    
    temp = data_1NET(1).temp; % just grab the first since all identical
    
    Sus = cat(1,data_1NET.Sus);
    %     Sus = mean(Sus,1);
    Sus = Sus./max(max(Sus(:,start:end)));
    
    Spec_Heat = cat(1,data_1NET.Spec_Heat);
    Spec_Heat = Spec_Heat./max(max(Spec_Heat(:,start:end)));

    
    Phi = cat(1,data_1NET.Phi);
    %     Phi = mean(Phi,1);
    Phi = Phi./max(max(Phi(:,start:end)));
    
    PhiSus = cat(1,data_1NET.PhiSus);
    PhiSus = PhiSus./max(PhiSus(:,start:end));
    
    Ener = cat(1,data_1NET.Ener);
    Ener = (Ener - min(Ener)) / ( max(Ener) - min(Ener) );
    % Ener = Ener./abs(max(Ener(:,start:end)));

    Mag = cat(1,data_1NET.Mag);
    
    figure(10)
    subplot(3,3,iNet)
    switch plotCase
        
        case 'phi'
            h1 = plot(temp,Sus,'.','Color',colors(1,:),'MarkerSize',17,'DisplayName','\chi');
            set(gca,'XScale','log')
            grid minor
            xticks(xt)
            hold all
            h2 = plot(temp,Phi,'+','Color',colors(2,:),'MarkerSize',mkSize,'DisplayName','\Phi', 'LineWidth',lWidth);
            set(gca,'XScale','log')
            grid minor
            xticks(xt)
            
            if drawPhiSus
                h6 = plot(temp,PhiSus,'s','Color',colors(8,:),'MarkerSize',mkSize,'DisplayName','\sigma^2_t(\Phi)','LineWidth',lWidth);
                set(gca,'XScale','log')
                grid minor
                xticks(xt)
            end
            
            if drawLine
                %%%%%%%%%%% Tmin Line %%%%%%%%%%%
                IndC = find(abs(temp - Out.Tmin_NET(iNet)) < 0.001);
                line([temp(IndC),temp(IndC)],[0,1],'LineStyle','--','Color','k','LineWidth',4)
                str = sprintf('$$T_{min} = %.2g $$',temp(IndC));
                text('Interpreter','latex','Position',[temp(IndC+3) 1/10],'String',str,'FontSize',fontSize)
                %%%%%%%%%%% Tc Line %%%%%%%%%%%
                if drawTcLine
                    IndC = find(abs(temp - Out.Tc_NET(iNet)) < 0.001);
                    line([temp(IndC),temp(IndC)],[0,Sus(IndC)],'Color','r')
                    str = sprintf('$$T_{c} = %.2g $$',temp(IndC));
                    text('Interpreter','latex','Position',[temp(IndC+1) 3/10],'String',str,'FontSize',fontSize)
                end
                if drawPhiSusLine
                    %%%%%%%%%%% Tcphi Line %%%%%%%%%%%
                    IndC = find(abs(temp - Out.Tcphi_NET(iNet)) < 0.001);
                    line([temp(IndC),temp(IndC)],[0,PhiSus(IndC)],'Color','g')
                    str = sprintf('$$T_{c,Phi} = %.2g $$',temp(IndC));
                    text('Interpreter','latex','Position',[temp(IndC+1) 5/10],'String',str,'FontSize',fontSize)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
            end
            
            if iNet == 3
                if drawPhiSus
                    lh = legend([h1,h2,h6]);
                    lh.FontSize = fontSize;
                else
                lh = legend([h1,h2]);
                lh.FontSize = fontSize;
                end
            end
            
        case 'ener'
            h3 = plot(temp,Ener,'.','Color',colors(3,:),'MarkerSize',mkSize,'DisplayName','E');
            hold all
            h4 = plot(temp,Mag,'.','Color',colors(4,:),'MarkerSize',mkSize,'DisplayName','M');
            
            if iNet == 1
                lh = legend([h3,h4]);                
                lh.FontSize = fontSize;
            end
            
        case 'sus/spec'
            h1 = plot(temp,Sus,'.','Color',colors(1,:),'MarkerSize',mkSize,'DisplayName','\chi');
            hold all
            h5 = plot(temp,Spec_Heat,'.','Color',colors(8,:),'MarkerSize',mkSize,'DisplayName','C');
            
            if iNet == 1
                lh = legend([h1,h5]);
                lh.FontSize = fontSize;
            end
    end
    
    
%     bl = boundedline(temp, mean(Sus,1), std(Sus),'.', ...
%         temp, mean(Phi,1), std(Phi),'.', ...
%         'cmap', colors);
    
    title([Networks{iNet}])
    
%     lh = legend(hl);
%     legnames = {'\chi', '\phi'};
%     for i = 1:length(legnames),
%         str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors(i, 1), colors(i, 2), colors(i, 3), legnames{i})];
%     end    
%     lh.String = str;
%     lh.FontSize = 15;
    
    % move a bit closer
%         lpos = lh.Position;
%         lpos(1) = lpos(1) + 0.15;
%         lh.Position = lpos;

    axis tight
    
    set(gca, ...
        'FontSize'    , 14        , ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.02 .02] , ...
        'XMinorTick'  , 'on'      , ...
        'YMinorTick'  , 'on'      , ...
        'YGrid'       , 'on'      , ...
        'XColor'      , [.3 .3 .3], ...
        'YColor'      , [.3 .3 .3], ...
        ...'YTick'       , 0:0.25:1.25, ...
        ...'YTickLabeL'  , [0:0.25:1.25],       ...
        'LineWidth'   , 1         );
   
    
%     if iNet == 6
%         ylim([0 1])
%     end
    
    
%     pause(0.1)
end
l1 = suplabel('Temperature','x',[.1 .10 .84 .84]);
l2 = suplabel('au','y',[.12 .1 .84 .84]);
l3 = suplabel('Brain Networks','t',[.1 .1 .84 .84]);

set(l1, 'FontSize',20)
set(l2, 'FontSize',20)
set(l3, 'FontSize',20)

% suptitle('\chi, \phi vs. Temperature')
% 
filename = 'figures_paper/mean_Phi_Sus_brain';
export_fig(filename, '-pdf', '-tiff', '-painters')

clearvars -except data_NET_met data_NET_glaub data_RAND_met data_RAND_glaub Out


%% Visualize Timeseries (NETS)

data = data_NET_met;

S = data(1).S;
Phi_Sub = data(1).Phi;
J = data(1).J;

corrWindow = 20;

for iTemp = 103:size(S,3)
    pause()
    Phi_T = Phi_Sub(iTemp);
    
    startFrame = 1;
    
    for iTime = 1:2:size(S,2)
        figure(1)
        subplot(2,2,1)
        imagesc(J)
        xticks([])
        yticks([])
        axis square
        title('Connectivity')
        
        subplot(2,2,2)
        bar(S(:,iTime,iTemp))
        axis square
        axis tight
        ylim([-2 2])
        xticks(1:size(S,1))
        yticks([])
        title(['\Phi = ',num2str(Phi_T) ,', t = ',num2str(iTime)])        
        
        subplot(2,2,3)
        
        if iTime > (startFrame + corrWindow)
            startFrame = startFrame + 1;
        end
        
        Corr_DTI = corrcoef(S(:,startFrame:iTime,iTemp)');
                
        imagesc(Corr_DTI)
        xticks([])
        yticks([])
        axis square
        title('Dynamic Correlations')
        
        subplot(2,2,4)
        
        Corr_DTI_full = corrcoef(S(:,1:iTime,iTemp)');
        
        imagesc(Corr_DTI_full)
        xticks([])
        yticks([])
        axis square
        title('Cumulative Correlations')
        
        pause(0.01)
    end
end

clearvars -except data_NET_met data_NET_glaub data_RAND_met data_RAND_glaub Out

%% Visualize the Random networks' correlations in order of the minimum distances first

homePhi = '../'; % Should lead to main 'Phi' directory
Corr_FMRI = load([homePhi, '/connectome_data/mean_struct_corr.mat']);
Corr_FMRI = Corr_FMRI.mean_corr_fc;
load([homePhi, 'reduced_networks/indices_84.mat'])

Networks = {'Aud', 'DMN', 'ECN_L', 'ECN_R', 'Salience', 'Sensorimotor', 'VISL', 'VISM', 'VISO'};

PD = Out.PD_RAND;

numExamples = 20;
numPerms = 120; % length(perms(1:5)) the number of permutations for a single Jij
permVecs = perms(1:5); % permutations matrix

for iNetwork = 1:length(Networks)
    
    Ind = indices.(Networks{iNetwork});
    
    dists = squeeze(PD(iNetwork,:,:));    
    
    % finding best Jij that minimizes distance
    [~, Isim] = sort(min(dists'));
    
    % find the min temps for each Jij
    
    % sort dists
    dists = dists(Isim,:);
    [~, Itemp] = min(dists');
    
    
    for iExample = 1:numExamples
        figure(1)
        subplot(1,3,1)
        % imagesc(Corr_FMRI(Ind,Ind))
        imagesc(atanh(Corr_FMRI(Ind,Ind)))
        axis  square
        caxis([0 1])
        title('\rho_{emp}')
        
        %find the correct Jij and rearrange appropriately
        simNum = Isim(iExample);
        JInd = floor(simNum/numPerms);
        
        permInd = mod(simNum,numPerms);
        if permInd == 0
            permInd = numPerms;
        end
        
        IJ = permVecs(permInd,:);
        J = data_RAND_met(JInd).J;
        J = J(IJ,IJ);
        
        Corr_DTI = data_RAND_met(JInd).Corr_DTI;
        Corr_DTI = Corr_DTI(IJ,IJ, Itemp(iExample));
        
        minDistr = dists(iExample,Itemp(iExample));
        
        subplot(1,3,2)
        
        % imagesc(Corr_DTI)
        imagesc(atanh(Corr_DTI))
        axis  square
        caxis([0 1])
        title('\rho_{sim,rand}')
        xlabel(['d = ', num2str(minDistr)])
        
%         subplot(1,3,3)
        
%         imagesc(J)
%         axis  square
%         caxis([-1 1])
%         title('J_{ij}')
                
        
        subplot(1,3,3)
        
        Corr_DTIn = data_NET_met(iNetwork).Corr_DTI;
        [minDistn, TInd] = min(Out.PD_NET(iNetwork,:,:)); 
        Corr_DTIn = squeeze(Corr_DTIn(:,:,TInd));
        
        % imagesc(Corr_DTIn)
        imagesc(atanh(Corr_DTIn))
        axis square
        caxis([0 1])
        title('\rho_{sim,brain}')
        xlabel(['d = ', num2str(minDistn)])
        
        titlestr = ['Sim# ', num2str(simNum),', J# ', num2str(JInd), ', Perm#: ',num2str(permInd)];
        
        suptitle(titlestr)
        
        suplabel(['Network = ', Networks{iNetwork}],'x')
        
        
        pause()
        
    end
    
end


clearvars -except data_NET_met data_NET_glaub data_RAND_met data_RAND_glaub Out

