%% Figure 2D
% Occurrences of BranchPair and Synapses vs. size of subvolume
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
figurePath = [matlabPath 'output\Figures\']; 
tablePath = [matlabPath 'output\Tables\']; 
filepath = [matlabPath 'data\branch_pairs_synapses\'];
            
%%
edgeLength = [100 50 25 10 5 1];
syn = []; % synapses
pairs = []; % branch pairs
g = []; 

meanPairs = nan(1,length(edgeLength));
meanSyn = nan(1,length(edgeLength));
sampleSize = nan(1,length(edgeLength));

for i = 1:length(edgeLength)
   cube_sz_str = num2str(edgeLength(i),'%d');
   filename = [filepath 'cube-stats_' cube_sz_str '-' ...
            cube_sz_str '-' cube_sz_str '.csv']; 
   tbl = readtable(filename);
   
   pairstmp = tbl.branches_pre.*tbl.branches_post;
   syntmp = tbl.boutons;
   idxDel = pairstmp==0; % Delete all cubes with no overlap
   pairstmp(idxDel) = [];
   syntmp(idxDel) = [];
   
   pairs = [pairs; pairstmp];
   syn = [syn; syntmp];    
   
   meanSyn(i) = mean(syntmp); 
   meanPairs(i) = mean(pairstmp);
   sampleSize(i) = sum(~idxDel);
   
   g = [g; edgeLength(i).*ones(size(syntmp))];
   
   fprintf('%d: %d (%d)\n',edgeLength(i),numel(idxDel),sum(idxDel));
end

%%
f1 = figure(1);
clf;
plot(edgeLength,meanPairs,'ko-','MarkerFaceColor','k','MarkerEdgeColor','none');
hold on;
plot(edgeLength,meanSyn,'ko-','MarkerFaceColor','k','MarkerEdgeColor','none');
xlabel('edge length (um)')
ylabel('occurrences');
set(gca,'Yscale','log','Xscale','log','Box','off','TickDir','out', ...
    'XTick',[1 10 100],'YTick',[10^0 10^3 10^6 10^9]);
% axis square;
set(f1,'PaperPositionMode', 'auto','Position',[0 0 300 200]);
print(f1,'-dsvg','-r600',[figurePath 'BranchPairs_vs_Synapses.svg']); 

% Difference in order magnitude
diffOrder = log10(meanPairs./meanSyn);
fprintf('Magnitude Order Difference: %.2f (%.2f %.2f)\n', ...
    median(diffOrder),min(diffOrder),max(diffOrder));

%% Table
tblResult = table;
tblResult.edgeLength = edgeLength';
tblResult.meanBranchPairs = meanPairs'; 
tblResult.meanSynapses = meanSyn';
tblResult.numVolumes = sampleSize'; 
writetable(tblResult,[tablePath 'PetersRule_BranchPairs_vs_Synapses.csv']); 