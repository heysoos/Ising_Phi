%% load connectivity matrix
wd = '/home/heysoos/Dropbox/Masters Thesis/Codes/Phi';
addpath('paper_figures/distinguishable_colors')
cd(wd);

J = load('MJ_all.mat');
J = J.MJ_all_mn;
J = J./(max(max(J)));

dist = distances(graph(1./J));
Jdist = 1./dist;
Jdist(isinf(Jdist)) = 0;
Jdist = Jdist./max(max(Jdist));

%% visualize graph
figure(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,3,1)

imagesc(J)
axis square
title('Brain Connectivity Map')
xlabel('ROIs')
ylabel('ROIs')
set(gca,'FontSize',18)
colormap(linspecer)
colorbar
xticklabels({})
yticklabels({})
set(gca,'xtick',[],'ytick',[])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,3,2)

Jthresh = J;
Jthresh(J < 0.006*max(max(J))) = 0;
G = graph(Jthresh);
deg = degree(G);

edgeC = [0.7, 0.7, 0.7];

LWidths = G.Edges.Weight/max(G.Edges.Weight);
LWidths = log(LWidths)-min(log(LWidths)) + eps; % normalized
H1 = plot(G, 'Layout', 'subspace', 'MarkerSize',log(deg+1+eps), ...               % node size in log scale
    ...'NodeCData',deg,...                             % node color by degree
    'EdgeColor',edgeC,'EdgeAlpha',0.3,'LineWidth',LWidths);
title('Projected Connectivity Graph')
set(gca,'FontSize',18)
axis square
xticklabels({})
yticklabels({})

H1.NodeLabel = {};
H1.NodeColor = 'black';

set(gca,'Color','w')
box off
axis off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,3,3)

imagesc(Jdist)
axis square
title('Inverse Distance Map')
xlabel('ROIs')
ylabel('ROIs')
set(gca,'FontSize',18)
colormap(linspecer)
colorbar
xticklabels({})
yticklabels({})
set(gca,'xtick',[],'ytick',[])                                
%% highlight networks, looped (zmaps, from pre-sorted data structure)

load('/home/heysoos/Dropbox/Masters Thesis/Codes/Phi/reduced_networks/reduced_Networks.mat')
% load('/home/heysoos/Dropbox/Masters Thesis/Codes/Phi/reduced_networks/indices.mat')

numNodesH = 15; % number of nodes to highlight
numNodesJ = 5; % number of nodes to extract for J_ij

% set up highlight colours
c = linspecer(length(reduced_Networks));

NET = 2 %1:9; % 1:9 for all networks, or pick a specific network

figure(3)
suptitle('Inverse Distance Map')
set(gca,'FontSize',20)

figure(4)
suptitle('Motifs')
set(gca,'FontSize',20)

figure(5)
suptitle('Connectivity Map')
set(gca,'FontSize',20)

for iNet = NET % highlight a particular network (see graphs_paper.m for all networks)
    
    
    % highlight each network's ROIs and edges    
    
    NETnodes = cat(1,...
        reduced_Networks(iNet).ROI_sort{end-(numNodesH-1):end,1}...
        );
    figure(1)
    subplot(1,3,2)
    title('Sub-Network Highlighted')
    for iROI = 1:length(NETnodes)
        
        NETedges = intersect(neighbors(...
            G,NETnodes(iROI)),NETnodes(NETnodes~=NETnodes(iROI)...
            ));
        highlight(H1, NETnodes(iROI),NETedges,'EdgeColor',c(iNet,:))
        highlight(H1, NETnodes(iROI),'NodeColor',c(iNet,:))
        
    end
    
    
    % generate subgraphs seperately for each network
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(3)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(3,3,iNet)
    
    % get the ROIs for each network from the connectivity distance
    Ind = cat(1,reduced_Networks(iNet).ROI_sort{end-4:end,1});
    
    imagesc(Jdist(Ind,Ind));
    axis square    
    txt1 = reduced_Networks(iNet).Network;
    title(txt1)
    xlabel('ROI')
    ylabel('ROIs')
    set(gca,'FontSize',18)
    colormap(linspecer)
    colorbar
    xticklabels({})
    yticklabels({})

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(3,3,iNet)
    
    % highlight each network's ROIs and edges  
    NETnodes = cat(1,reduced_Networks(iNet).ROI_sort{end-(numNodesJ-1):end,1});
    
    Hsub = subgraph(G,NETnodes);
    Hsub1 = plot(Hsub,'Layout','force');    
    txt1 = reduced_Networks(iNet).Network;
    title(txt1)
    axis square
    Hsub1.NodeLabel = {};
    Hsub1.NodeColor = 'black';
    set(gca,'FontSize',18)
    xticklabels({})
    yticklabels({})
    axis square
    axis off
    box off
    set(gca,'xtick',[],'ytick',[])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(5)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(3,3,iNet)
    
    % get the ROIs for each network from the connectivity distance
    Ind = cat(1,reduced_Networks(iNet).ROI_sort{end-4:end,1});
    
    imagesc(J(Ind,Ind));
    axis square
    
    txt1 = reduced_Networks(iNet).Network;
    title(txt1)    
    xlabel('ROI')
    ylabel('ROIs')
    set(gca,'FontSize',18)
    colormap(linspecer)
    colorbar
    xticklabels({})
    yticklabels({})
    
    pause(0.01)
end