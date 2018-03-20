%% load connectivity matrix
wd = '/home/heysoos/Dropbox/Masters Thesis/Codes/Phi';
addpath('paper_figures/distinguishable_colors')
cd(wd);

J = load('MJ_all.mat');
J = J.MJ_all_mn;
J = J./(max(max(J)));

%% visualize graph

Jthresh = J;
Jthresh(J < 0.006*max(max(J))) = 0;
G = graph(Jthresh);
deg = degree(G);

edgeC = [0.7, 0.7, 0.7];

LWidths = G.Edges.Weight/max(G.Edges.Weight);
LWidths = log(LWidths)-min(log(LWidths)) + eps; % normalized

subplot(1,2,1)


imagesc(J)
axis square
title('J_{ij} Whole Brain Connectivity')
xlabel('ROIs')
colormap(linspecer)

subplot(1,2,2)

% colormap cool                                       % set color map
H1 = plot(G, 'Layout', 'subspace', 'MarkerSize',log(deg+1+eps), ...               % node size in log scale
    ...'NodeCData',deg,...                             % node color by degree
    'EdgeColor',edgeC,'EdgeAlpha',0.3,'LineWidth',LWidths);
title('J_{ij }Whole Brain Graph with Highlighted Networks')

H1.NodeLabel = {};
H1.NodeColor = 'black';

axis square



%% highlight networks (zmaps, from pre-sorted data structure)

load('/home/heysoos/Dropbox/Masters Thesis/Codes/Phi/reduced_networks/reduced_Networks.mat')

numNodes = 10; % number of nodes assigned to each network

% set up highlight colours
c = linspecer(length(reduced_Networks));

% highlight each network's ROIs and edges
for iNet = 2%:length(reduced_Networks)
    
    NETnodes = cat(1,...
        reduced_Networks(iNet).ROI_sort{end-(numNodes-1):end,1}...
        );
    
    for iROI = 1:length(NETnodes)
        
        NETedges = intersect(neighbors(...
            G,NETnodes(iROI)),NETnodes(NETnodes~=NETnodes(iROI)...
            ));
        highlight(H1, NETnodes(iROI),NETedges,'EdgeColor',c(iNet,:))
        highlight(H1, NETnodes(iROI),'NodeColor',c(iNet,:))
        
    end
    %     pause()
end

%% generate subgraphs seperately for each network

load('/home/heysoos/Dropbox/Masters Thesis/Codes/Phi/reduced_networks/reduced_Networks.mat')

numNodes = 10; % number of nodes assigned to each network

% set up highlight colours
c = linspecer(length(reduced_Networks));

% highlight each network's ROIs and edges
for iNet = 1:length(reduced_Networks)
    
    NETnodes = cat(1,reduced_Networks(iNet).ROI_sort{end-(numNodes-1):end,1});
    
    Hsub = subgraph(G,NETnodes);
    plot(Hsub,'Layout','subspace')
    pause()
    
end