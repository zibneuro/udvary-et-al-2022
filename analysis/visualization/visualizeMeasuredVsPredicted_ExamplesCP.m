%% Figure 6C
% Histograms of example connection probably distributions vs observed ones
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

%%
matlabPath = 'D:\udvary-et-al-2022\analysis\';
figurePath = [matlabPath 'output\Figures\']; 
filePath = [matlabPath 'data\cellular_connectivity\stats\'];

%%
tbl = readtable([filePath 'summary.csv']); 
CP_measured = tbl.empirical'.*100;
CP_AVG_predicted = tbl.avg_slices_merged'.*100;
CP_SD_predicted = tbl.std_slices_merged'.*100;

% Examples for paper
% 3: VPM -> L4 Bruno & Sakman 2006 (red)
% 46: L4sp->L4sp Sun et al. ,2006 (anatomy EM) (green)
% 36: L4->L23 (Yoshimura, Optical Stimulation)  (blue)
IDList = [3 46 36];
col = 'rgbc'; 
xBins = [0:1:100]; % 1% bins

f1 = figure(1);
clf; 

for i = 1:length(IDList)
    prefix = tbl.descriptor(IDList(i));
    filename = [filePath prefix{:} '_histogram.txt'];
    tbltmp = readtable(filename);
    y = tbltmp.Var1'; 
    y = y./max(y); 
    subplot(3,1,i);
%     bar(xBins,y,1,'EdgeColor','none','FaceColor','k');
    stairs([xBins(1) xBins(1:end) xBins(end)],[0 y 0], ...
                'Color','k'); 
    hold on;
    plot(CP_measured(IDList(i)),1,'o','MarkerEdgeColor','none', ...
        'MarkerFaceColor','g'); 
    plot(CP_AVG_predicted(IDList(i)),1,'o','MarkerEdgeColor','none', ...
        'MarkerFaceColor','k'); 

    % Range of M +- SD
    rangeSD = CP_AVG_predicted(IDList(i)) + ...
                CP_SD_predicted(IDList(i)).*[1 -1];
    rangeSD(rangeSD<0) = 0;
    rangeSD(rangeSD>100) = 100;
    
    % Mean and range of SD
    plot(rangeSD,[1 1],'k-','LineWidth',1);
    for j = 1:2
        plot([rangeSD(j) rangeSD(j)],1+[-0.05 0.05],'k-'); 
    end
   
   set(gca,'Box','off','TickDir','out','XLim',[-0.5 100], ...
       'XTick',[0:25:100]);
end

%%
set(f1,'PaperPositionMode', 'auto','Position',[0 0 200 300]);
print(f1,'-dsvg','-r600',[figurePath 'MeasuredVsPredicted_ExampleCP.svg']);