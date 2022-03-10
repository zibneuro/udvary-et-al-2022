%% Figure 3G
% Visualize Deviation of Fully Recurrent and Chain motifs
% Deviation p/p_rand vs. neurons per motif
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
figurePath = [matlabPath 'output\Figures\']; 
tablePath = [matlabPath 'output\Tables\']; 
postfix = 'C2'; 

%% Chain Motifs
numTrials = 10000;
numNodesListChains = [3 4 5 6 7 8 9 10]; 
devChains = nan(numel(numNodesListChains),1); 
xStr = cell(1,numel(numNodesListChains)); 
p1 = nan(size(devChains)); 
p2 = nan(size(devChains)); 

for numNodes = numNodesListChains
    
    load([matlabPath 'data\higherOrder\' ...
            'ChainMotif_NumTrials_' num2str(numTrials) ...
            '_NumCells_' num2str(numNodes) '_' postfix '.mat']);
    
    p = mean(edgeStats.pEdgeMotif,2)';
    p1(numNodes==numNodesListChains) = mean(p);
    p2(numNodes==numNodesListChains) = edgeStats.pEdgeUniformMotif;
    devChains(numNodes==numNodesListChains) = mean(p)./edgeStats.pEdgeUniformMotif;
    xStr{numNodes==numNodesListChains} = num2str(numNodes); 
end

t = table(numNodesListChains',p1,p2,devChains, ...
            'VariableNames',{'NumNodes','p(model)','p(random)','deviation'});
writetable(t,[tablePath 'ChainMotifs_' postfix '.csv']);

%% Fully Recurrent Motifs
numTrials = 10000000;
load([matlabPath 'data\higherOrder\' ...
    'FullyRecurrentMotif_NumTrials_' num2str(numTrials) '_' postfix '.mat']);
devFullyRec = nan(numel(numNodesListFullyRec),1); 

for numNodes = numNodesListFullyRec
    
    
    devFullyRec(numNodes==numNodesListFullyRec) = ...
        pRec(numNodes==numNodesListFullyRec)./pRecRandom(numNodes==numNodesListFullyRec);
    
    fprintf('%d: %.2e vs %.2e = %.2e\n',numNodes, ...
            pRec(numNodes==numNodesListFullyRec), ...
            pRecRandom(numNodes==numNodesListFullyRec), ...
            devFullyRec(numNodes==numNodesListFullyRec));

end

t = table(numNodesListFullyRec',pRec',pRecRandom',devFullyRec, ...
            'VariableNames',{'NumNodes','p(model)','p(random)','deviation'});
writetable(t,[tablePath 'RecurrentMotifs_' postfix '.csv']);

%% Plot
if sum(numNodesListChains==numNodesListFullyRec)~=numel(numNodesListFullyRec)
    error('Not same list of number of nodes');
end

f1 = figure(1);
clf;
for i = 1:2
    subplot(2,1,i)
    xLimPlot = [0 numel(numNodesListChains)]+0.5; 
%     plot(xLimPlot,[1 1],'-','Color',[0.7 0.7 0.7]); 
%     hold on; 
    if i==1
        plot(devFullyRec,'k-o','MarkerFaceColor','k'); 
        yLimRange = [1 1e26];
        yTick = [1 1e13 1e26];
    else
        plot(devChains,'k-o','MarkerFaceColor','k'); 
        yLimRange = [0.15 1];
        yTick = [0.25 0.5 1];
    end
    xlim(xLimPlot);
    set(gca,'Box','off','TickDir','out','XTick',[1:numel(numNodesListChains)], ...
        'XTickLabel',xStr,'YScale','log','YLim',yLimRange,'YTick',yTick);
%     set(gca,'Box','off','TickDir','out','XTick',[1:numel(numNodesListChains)], ...
%         'XTickLabel',xStr,'YScale','log');
    ylabel('deviation');
    xlabel('#nodes');
end
set(f1,'PaperPositionMode', 'auto','Position',[0 0 350 300]);
print(f1,'-dsvg','-r600',[figurePath 'Deviation_Chains_Rec_NumNodes_Log.svg']);