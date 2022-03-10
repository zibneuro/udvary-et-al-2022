%% Figure 2H
% Occurrences of neuron pairs vs. connected overlap volumes
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
addpath(genpath([matlabPath 'functions\']));
figurePath = [matlabPath 'output\Figures\']; 
tablePath = [matlabPath 'output\Tables\']; 
storagePath = [matlabPath 'data\dataPetersRuleBatch\'];

%%
batchID_list = [1 2 3 4 5 6 7 8 9];

load([storagePath 'dsc_overlapping-cubes_batch-0.mat']);
Npairs = result.Npairs;
edge_sz = result.edge_sz;
perc_edges = result.perc_edges;
perc_p_mean = result.perc_p_mean;
perc_p_sd = result.perc_p_sd;
perc_p_N = result.perc_p_N;
numBins = length(perc_p_mean);

for batchID = batchID_list
    load([storagePath 'dsc_overlapping-cubes_batch-' num2str(batchID) ...
                '.mat']);
    Npairs = Npairs + result.Npairs;

    for j = 1:numBins
        [M,SD] = totalMSD([perc_p_mean(j) result.perc_p_mean(j)], ...
                            [perc_p_sd(j) result.perc_p_sd(j)], ...
                            [perc_p_N(j) result.perc_p_N(j)]);
        perc_p_sd(j) = SD;
        perc_p_mean(j) = M;
        perc_p_N(j) = perc_p_N(j) + result.perc_p_N(j);
    end
end

%%
% create dedicated 0 and 1 bin
x = [0 (perc_edges(1:end-1)+edge_sz/2) 1].*100;
y_mean = perc_p_mean.*perc_p_N; 

% Calculate SEM
y_sem = perc_p_sd./sqrt(perc_p_N);
y_sem = y_sem.*perc_p_N;
y_min = y_mean-y_sem;
y_max = y_mean+y_sem;

fprintf('#pairs = %.0f\n',Npairs);
fprintf('  0%%: %.0f \n',y_mean(1));
fprintf('  100%%: %.0f \n',y_mean(end));

%%
smoothingFac = 0.9;
f1 = figure(1);
clf;
subplot(1,2,1);
y_mean_smooth = smoothdata(y_mean,'movmedian','SmoothingFactor',smoothingFac);
hold on;
y_min_smooth = smoothdata(y_min,'movmedian','SmoothingFactor',smoothingFac);
y_max_smooth = smoothdata(y_max,'movmedian','SmoothingFactor',smoothingFac);
plot(x,y_mean_smooth,'k-','LineWidth',1);
set(gca,'YScale','log','TickDir','out','Box','off','YLim',[0.5 Npairs], ...
    'XTick',[0:25:100],'YTick',[1 1e2 1e4 1e6 1e8]);
ylabel('occurrences');
xlabel('overlap %');
title('smoothed');

subplot(1,2,2);
plot(x,y_mean,'k-','LineWidth',1);
hold on;
set(gca,'YScale','log','TickDir','out','Box','off','YLim',[0.9 Npairs], ...
    'XTick',[0:25:100],'YTick',[1 1e2 1e4 1e6 1e8]);
ylabel('occurrences');
xlabel('overlap %');
title('raw');

%%
tblResult = table;
tblResult.ConnectionsPerCellPair_in_overlap_percentage = x';
tblResult.MeanOccurrencesSmoothed = y_mean_smooth';
tblResult.MinOccurrencesSmoothed = y_min_smooth';
tblResult.MaxOccurrencesSmoothed = y_max_smooth';
tblResult.MeanOccurrences = y_mean';
tblResult.MinOccurrences = y_min';
tblResult.MaxOccurrences = y_max';
writetable(tblResult,[tablePath 'OccurrencesPerOverlap.csv']); 

%%
set(f1,'PaperPositionMode','auto','Position',[0 0 600 200]); 
print(f1,'-painters','-dsvg','-r600',[figurePath ...
                    'PetersRule_Conn_Overlap.svg']); 
saveas(f1,[figurePath 'PetersRule_Conn_Overlap.jpg']);