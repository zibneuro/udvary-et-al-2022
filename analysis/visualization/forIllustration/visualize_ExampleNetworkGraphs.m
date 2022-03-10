%% Figure 3D/C
% Extract example network of 50 nodes that has the 3 example cells in it
% only for C2 column
% only for illustration purposes!
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
figuresPath = [matlabPath 'output\Figures\']; 
storagePath = [matlabPath 'data\'];
tablePath = [matlabPath 'output\Tables\'];
addpath(genpath([matlabPath 'functions\']));
motifsPath = [storagePath 'motifs.csv'];

%% When importing in illustrator:
% Delete white background nodes!
% Delete "stroke" and do proper stroke!
% Arrows are weirdly cropped and isolated ... 

%%
L5PT_ID = 301854;
L2PY_ID = 748854;
L6ACC_ID = 199678;
load([storagePath 'NetworkGraphExample.mat'], ...
            'pMatrix','exampleIDs','exampleCTs');
diagID = logical(eye(size(pMatrix))); 
p_avg = mean(pMatrix(~diagID)); 

% Plot Settings
% Color Settings
col = [0 0.71 0.14; 0.86 0.42 0.09; 0.01 0.01 0.01];  
% Cell Type colors L5PT L2PY L6ACC
grayCol = [0.5 0.5 0.5]; 
nodeColHidden = [1 1 1]; 
edgeHighCol = 'k'; 
% nodeHighCol = 'r'; 
% nodeStrokeCol = 'k'; 
rndGraphCol = 'r'; 
nodeColRnd = [1 0.7 0.7]; 

% Size settings
mrksz = 11; 
arrsz = 10; 
edgeHighWdth = 1; 
lineWdth = 0.25; 
arrwPos = 0.99; 

% Hardcoded example subgraph
nodeHighlighted = [find(exampleIDs==L5PT_ID) ...
                    find(exampleIDs==L2PY_ID) ...
                    find(exampleIDs==L6ACC_ID)];
exampleMotifID = 10;
exampleMotifNumEdges = 3; 

%%
rng(2984924029);

% Each Network should contain numEdges (edges based on average p)
numEdges = round((numel(pMatrix)-size(pMatrix,1))*p_avg);

% Get realization
ptmp = sort(pMatrix(:),'descend');
m1 = pMatrix>ptmp(numEdges+1);

% Shuffle m1 (same as taking p_avg)
valtmp = m1(~diagID);
idxRnd = randperm(numel(valtmp)); 
mRand = false(size(m1));
mRand(~diagID) = valtmp(idxRnd);

% Generate network 2 with the same number of edges
% basically generate random networks until numEdges is matched
m2 = false(size(m1));
while sum(m2(:))~=numEdges
    r = rand(size(pMatrix)); 
    m2 = pMatrix>r; 
end

% Calculate probability of each realization (in log space)
% p1 = prod(pMatrix(m1 & ~diagID))*prod(1-pMatrix(~m1 & ~diagID));
% p2 = prod(pMatrix(m2 & ~diagID))*prod(1-pMatrix(~m2 & ~diagID));
p1 = sum(log10(pMatrix(m1 & ~diagID))) + sum(log10(1-pMatrix(~m1 & ~diagID)));
p2 = sum(log10(pMatrix(m2 & ~diagID))) + sum(log10(1-pMatrix(~m2 & ~diagID)));
pRnd = sum(log10(pMatrix(mRand & ~diagID))) + sum(log10(1-pMatrix(~mRand & ~diagID)));
mImpossibleEdges = (mRand>0 & pMatrix==0);
mPossibleEdges = (mRand>0 & pMatrix>0);

if p1>p2
    fprintf('Network 1 more likely than Network 2\n');
else
    fprintf('Network 2 more likely than Network 1\n');
end
fprintf(' >> log10(p(Network1)) = %.4f\n >> log10(p(Network2)) = %.4f\n',p1,p2);
fprintf(' >> log10(p(NetworkRandom)) = %.4f\n',pRnd);
fprintf('Number of edges:\n >> Network 1: %d\n >> Network 2: %d\n >> NetworkRnd: %d\n', ...
        sum(m1(:)),sum(m2(:)),sum(mRand(:)));
fprintf(' >> Number of structurally impossible edges in random Network: %d\n', ...
        sum(mImpossibleEdges(:)));

% Save as matrix
writeNetworkAsCSV([tablePath 'exampleNetworks.csv'],pMatrix,m1,m2,mRand,exampleIDs); 
    
%% Check
% Go through all triplet motifs and find exemplary motif   
[exampleNodes1, exampleNodesRnd] = getNodesThatMatchMotif(m1,mRand,exampleMotifID,motifsPath);
[exampleNodes2] = getNodesThatMatchMotif(m2,mRand,exampleMotifID,motifsPath);
tmp = sort(exampleNodes1,2);
n1 = size(unique(tmp,'rows'),1);
tmp = sort(exampleNodes2,2);
n2 = size(unique(tmp,'rows'),1);
tmp = sort(exampleNodesRnd,2);
nRnd = size(unique(tmp,'rows'),1);
fprintf(['Motif ID %d:\n >> %d in Network1\n >> %d in Network2\n' ...
        ' >> %d in Random Network\n'], exampleMotifID,n1,n2,nRnd);

%%
% Create Digraph
G1 = digraph(m1); 
G2 = digraph(m2); 
Grnd = digraph(mRand);
Grnd_Imp = digraph(mImpossibleEdges);
Grnd_Po = digraph(mPossibleEdges);

%% Plotting
f1 = figure(1);
clf;

%% Barrel Cortex Graph 1
subplot(1,3,1); 
load([storagePath 'NodePositions.mat'],'X','Y');
% pplot1 = plot(G1,'Layout','force','NodeLabel','','NodeColor',nodeColHidden, ...
%                 'EdgeColor',grayCol,'EdgeAlpha',1,'LineWidth',lineWdth, ...
%                 'ArrowPosition',arrwPos,'ArrowSize',arrsz,'MarkerSize',mrksz); 
pplot1 = plot(G1,'XData',X,'YData',Y,'NodeLabel','','NodeColor',nodeColHidden, ...
                'EdgeColor',grayCol,'EdgeAlpha',1,'LineWidth',lineWdth, ...
                'ArrowPosition',arrwPos,'ArrowSize',arrsz,'MarkerSize',mrksz); 
hold on; 

% Highlight example nodes
for i = 1:length(X)         
    faceColtmp = grayCol.*1.5; 
    idxtmp = find(i==nodeHighlighted);
    
    if ~isempty(idxtmp)
        faceColtmp = col(idxtmp,:);
    end
    plot(X(i),Y(i),'o','MarkerEdgeColor','none', ...
            'MarkerFaceColor',faceColtmp,'MarkerSize',mrksz);
end

% Highlight example motifs
highlightInGraphEdgeBased(G1,exampleNodes1(1,:),1:length(X),pplot1, ...
                               edgeHighCol,edgeHighWdth);
xl = xlim;
yl = ylim; 
set(gca,'XTick',[],'YTick',[],'XLim',xl,'Ylim',yl);
axis square;

%% Barrel Cortex Graph 2
subplot(1,3,2); 
pplot2 = plot(G2,'XData',X,'YData',Y,'NodeLabel','','NodeColor',nodeColHidden, ...
                'EdgeColor',grayCol,'EdgeAlpha',1,'LineWidth',lineWdth, ...
                'ArrowPosition',arrwPos,'ArrowSize',arrsz,'MarkerSize',mrksz);
hold on; 

% Highlight example nodes
for i = 1:length(X)         
    faceColtmp = grayCol.*1.5; 
    idxtmp = find(i==nodeHighlighted);
    
    if ~isempty(idxtmp)
        faceColtmp = col(idxtmp,:);
    end
    plot(X(i),Y(i),'o','MarkerEdgeColor','none', ...
            'MarkerFaceColor',faceColtmp,'MarkerSize',mrksz);
end

% Highlight example motifs
highlightInGraphEdgeBased(G2,exampleNodes2(1,:),1:length(X),pplot2, ...
                                edgeHighCol,edgeHighWdth);
set(gca,'XTick',[],'YTick',[],'XLim',xl,'Ylim',yl);
axis square;

%% Random Graph
subplot(1,3,3); 
pplotrnd = plot(Grnd,'XData',X,'YData',Y,'NodeLabel','','NodeColor',nodeColHidden, ...
                'EdgeColor',rndGraphCol,'EdgeAlpha',1,'LineWidth',lineWdth, ...
                'ArrowPosition',arrwPos,'ArrowSize',arrsz,'MarkerSize',mrksz);
hold on; 

% Highlight example nodes
for i = 1:length(X)         
    faceColtmp = nodeColRnd; 
    idxtmp = find(i==nodeHighlighted);
    
    if ~isempty(idxtmp)
        faceColtmp = col(idxtmp,:);
    end
    plot(X(i),Y(i),'o','MarkerEdgeColor','none', ...
            'MarkerFaceColor',faceColtmp,'MarkerSize',mrksz);
end

% Highlight example motifs
highlightInGraphEdgeBased(Grnd,exampleNodesRnd(1,:),1:length(X),pplotrnd, ...
                                edgeHighCol,edgeHighWdth);
set(gca,'XTick',[],'YTick',[],'XLim',xl,'Ylim',yl);
axis square;

%%
set(f1,'PaperPositionMode','auto','Position',[0 0 1600 1000]); 
print(f1,'-painters','-dsvg','-r600',[figuresPath 'ExampleNetworkGraph_v2.svg']); 

%% Impossible vs Possible Edges
f2 = figure(2);
clf;

subplot(1,3,2); 
pplotrnd = plot(Grnd,'XData',X,'YData',Y,'NodeLabel','','NodeColor',nodeColHidden, ...
                'EdgeColor',rndGraphCol,'EdgeAlpha',1,'LineWidth',lineWdth, ...
                'ArrowPosition',arrwPos,'ArrowSize',arrsz,'MarkerSize',mrksz);
hold on; 

% Highlight example nodes
for i = 1:length(X)         
    faceColtmp = nodeColRnd; 
    idxtmp = find(i==nodeHighlighted);
    
    if ~isempty(idxtmp)
        faceColtmp = col(idxtmp,:);
    end
    plot(X(i),Y(i),'o','MarkerEdgeColor','none', ...
            'MarkerFaceColor',faceColtmp,'MarkerSize',mrksz);
end

% Highlight example motifs
highlightInGraphEdgeBased(Grnd,exampleNodesRnd(1,:),1:length(X),pplotrnd, ...
                                edgeHighCol,edgeHighWdth);
set(gca,'XTick',[],'YTick',[],'XLim',xl,'Ylim',yl);
axis square;

subplot(1,3,3); 
pplotrnd = plot(Grnd_Imp,'XData',X,'YData',Y,'NodeLabel','','NodeColor',nodeColHidden, ...
                'EdgeColor',rndGraphCol,'EdgeAlpha',1,'LineWidth',lineWdth, ...
                'ArrowPosition',arrwPos,'ArrowSize',arrsz,'MarkerSize',mrksz);
set(gca,'XTick',[],'YTick',[],'XLim',xl,'Ylim',yl);
axis square;
set(f2,'PaperPositionMode','auto','Position',[0 0 1600 1000]); 
print(f2,'-painters','-dsvg','-r600',[figuresPath 'ExampleNetworkGraph_v2_ImpossibleEdges.svg']); 

%% Functions
% Extract Nodes that are part of target motif
function [nodeList, nodeListRnd] = getNodesThatMatchMotif(m,mRand,motifID,motifsPath)
    
    % Get Motifs
    MotifSpectrum = getMotifSpectrum(motifsPath); 
    newIdxOrder = [1 3 5 6 4 7 2 8 10 11 12 13 9 14 15 16];
    idxMotif = find(MotifSpectrum.Class==newIdxOrder(motifID))';
    
    nodeList = []; 
    nodeListRnd = []; 
    sz = size(m,1); 
    
    for x = 1:sz
        for y = 1:sz
            
            if x==y
                continue;
            end
            
            for z = 1:sz
                
                if z==y || z==x
                    continue;
                end
                
                % Get Triplet motif
                triplet = zeros(3,3); 
                triplet(1,2:3) = m(x,[y z]); 
                triplet(2,[1 3]) = m(y,[x z]);
                triplet(3,1:2) = m(z,[x y]);
                
                tripletRand = zeros(3,3); 
                tripletRand(1,2:3) = mRand(x,[y z]); 
                tripletRand(2,[1 3]) = mRand(y,[x z]);
                tripletRand(3,1:2) = mRand(z,[x y]);
                
                for mm = idxMotif
                    tripletTarget = ...
                            squeeze(MotifSpectrum.Connectivity(mm,:,:));
                    if sum(triplet(:)==tripletTarget(:))==9
                        nodeList = [nodeList; x y z];
                    end
                    if sum(tripletRand(:)==tripletTarget(:))==9
                        nodeListRnd = [nodeListRnd; x y z];
                    end 
                end
            end
        end
    end
end

function [edgesExample] = highlightInGraphNodeBased(G,nodes,pplot,edgeColor,lineWdth)

    edgesExample = []; 
    for i = 1:length(nodes)
        for j = 1:length(nodes)
            if i==j
                continue;
            end
            idxOut = findedge(G,nodes(i),nodes(j)); 
            if idxOut~=0
                edgesExample = [edgesExample idxOut];
            end
        end
    end

    highlight(pplot,G.Edges.EndNodes(edgesExample,1), ...
                    G.Edges.EndNodes(edgesExample,2), ...
                    'EdgeColor',edgeColor,'LineWidth',lineWdth); 
end

function [edgesExample] = highlightInGraphEdgeBased(G,nodeList, ...
                                nodeHighlighted,pplot,edgeColor,lineWdth)
    
    edgesExample = []; 
    for ii = 1:size(nodeList,1)
        idxtmp = ismember(nodeHighlighted,nodeList(ii,:)); 
        
        if sum(idxtmp)==3
            nodestmp = find(idxtmp); 
            for i = 1:length(nodestmp)
                for j = 1:length(nodestmp)
                    if i==j
                        continue;
                    end
                    idxOut = findedge(G,nodestmp(i),nodestmp(j)); 
                    if idxOut~=0
                        edgesExample = [edgesExample idxOut];
                    end
                end
            end
        end
    end
    
    highlight(pplot,G.Edges.EndNodes(edgesExample,1), ...
                    G.Edges.EndNodes(edgesExample,2), ...
                    'EdgeColor',edgeColor,'LineWidth',lineWdth);   
end

function [] = writeNetworkAsCSV(filename,pMatrix,m1,m2,mRand,exampleIDs)
    fid = fopen(filename,'w');
    if fid==-1
       error('Could not open %s!',filename); 
    end
    
    fprintf(fid,['preID,postID,connection_probability,' ...
                    'network1,network2,random_network,\n']);
    for i = 1:numel(exampleIDs)
        for j = 1:numel(exampleIDs)
            fprintf(fid,'%d,%d,%.4f,%d,%d,%d,\n', ...
                    exampleIDs(i),exampleIDs(j), ...
                    pMatrix(i,j),m1(i,j),m2(i,j),mRand(i,j)); 
        end
    end
    fclose(fid);
end