%% Figure 3E
% Triplet Motif Deviatio
% - all cells
% - L5EXC cells in C2
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
figurePath = [matlabPath 'output\Figures\']; 
tablePath = [matlabPath 'output\Tables\']; 
load([matlabPath 'data\TripletMotifs_All.mat'],'tbl'); 

%% All vS1
idx = strcmp(tbl.A,'ALL-ALL') & strcmp(tbl.B,'ALL-ALL') & ...
            strcmp(tbl.C,'ALL-ALL');
p1Motif = tbl.pMotif(idx,:);
p1MotifRnd = tbl.pMotifRnd(idx,:);
dev1 = p1Motif./p1MotifRnd;
dev1log = log10(dev1);
dev1logNorm = dev1log./max(dev1log); 

%% Layer Example
idx = strcmp(tbl.A,'C2-L5EXC') & strcmp(tbl.B,'C2-L5EXC') & ...
            strcmp(tbl.C,'C2-L5EXC');
p2Motif = tbl.pMotif(idx,:);
p2MotifRnd = tbl.pMotifRnd(idx,:);
dev2 = p2Motif./p2MotifRnd;
dev2log = log10(dev2); 
dev2logNorm = dev2log./max(dev2log); 

%% Figure
f2 = figure(2);
clf; 
bar(1:15,[dev1logNorm(1:15)' dev2logNorm(1:15)'],1,'BaseValue',0,'EdgeColor','none'); 
set(gca,'TickDir','out','Box','off','XLim',[0.5 15.5],'XTick', ...
        [1:16],'YScale','linear'); 
set(f2,'PaperPositionMode', 'auto','Position',[0 0 350 200]);
print(f2,'-dsvg','-r600',[figurePath 'TripletMotifs_Examples.svg']);

%%
tblResult = table;
tblResult.motifID = [1:15]';
tblResult.probability_vS1 = p1Motif(1:15)'; 
tblResult.probability_random_vS1 = p1MotifRnd(1:15)'; 
tblResult.deviation_vS1 = dev1(1:15)'; 
tblResult.deviation_vS1_log = dev1log(1:15)'; 
tblResult.deviation_vS1_log_maxnormalized = dev1logNorm(1:15)';
tblResult.probability_L5_C2 = p2Motif(1:15)'; 
tblResult.probability_random_L5_C2 = p2MotifRnd(1:15)'; 
tblResult.deviation_L5_C2 = dev2(1:15)'; 
tblResult.deviation_L5_C2_log = dev2log(1:15)'; 
tblResult.deviation_L5_C2_log_maxnormalized = dev2logNorm(1:15)'; 

writetable(tblResult,[tablePath 'TripletExamples_vS1_and_L5_C2.csv']);