%% Figure 2H
% Overlapping and connected cell pairs vs. size of subvolume
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
figurePath = [matlabPath 'output\Figures\']; 
tablePath = [matlabPath 'output\Tables\']; 
filepath = [matlabPath 'data\overlapping_connected\'];

%%
edgeLength = [100 50 25 10 5 1];
meanOverlapPairs = nan(1,length(edgeLength));
meanConnPairs = nan(1,length(edgeLength));

for i = 1:length(edgeLength)
   cube_sz_str = num2str(edgeLength(i),'%d');
   filename = [filepath 'overlapping_connected_' cube_sz_str '-' ...
                    cube_sz_str '-' cube_sz_str '.csv']; 
   tbl = readtable(filename);
      
   meanConnPairs(i) = tbl.connected_pairs; 
   meanOverlapPairs(i) = tbl.overlapping_pairs;   
   
   fprintf('%d,%.0f,\n',edgeLength(i),tbl.overlapping_pairs);
end

%%
f1 = figure(1);
clf;
plot(edgeLength,meanConnPairs,'ko-','MarkerFaceColor', ...
            'k','MarkerEdgeColor','none');
hold on;
plot(edgeLength,meanOverlapPairs,'ko-','MarkerFaceColor', ...
            'k','MarkerEdgeColor','none');
xlabel('edge length (um)')
ylabel('occurrences');
set(gca,'Yscale','log','Xscale','log','Box','off','TickDir','out', ...
    'XTick',[1 10 100],'YTick',[1e7 1e8 1e9],'XLim',[1 100], ...
    'YLim',[1e7 1e9]);
set(f1,'PaperPositionMode', 'auto','Position',[0 0 300 200]);
print(f1,'-dsvg','-r600',[figurePath 'Overlap_vs_Connected.svg']); 

%% Table
tblResult = table;
tblResult.edgeLength = edgeLength';
tblResult.meanConnectedCellPairs = meanConnPairs'; 
tblResult.meanOverlappingCellPairs = meanOverlapPairs';
writetable(tblResult,[tablePath 'PetersRule_Overlap_vs_ConnectedCellPairs.csv']); 