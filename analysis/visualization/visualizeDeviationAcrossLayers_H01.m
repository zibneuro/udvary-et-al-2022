% Figure 7C
% Predicted vs. observed motif deviations in H01
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
figurePath = [matlabPath 'output\Figures\']; 
prefix = 'ALL';
datapath = [matlabPath 'data\H01\motifs_' prefix '\h01-subpopulations\'];
tblMotifs = readtable([datapath 'probabilities_16_h01-subpopulations.csv']);
sample = 200;
layers = 2:6; 
colGreen = [20 205 20]./255; 
colLayer = [107 62 152; 246 153 153; 27 175 72; 250 162 45; 31 120 180]./255;

%%
tbl = table();
tbl.A = tblMotifs.A;
tbl.B = tblMotifs.B;
tbl.C = tblMotifs.C; 
tbl.pMotif = [tblMotifs.motif_1_model tblMotifs.motif_2_model ...
              tblMotifs.motif_3_model tblMotifs.motif_4_model ...
              tblMotifs.motif_5_model tblMotifs.motif_6_model ...
              tblMotifs.motif_7_model tblMotifs.motif_8_model ...
              tblMotifs.motif_9_model tblMotifs.motif_10_model ...
              tblMotifs.motif_11_model tblMotifs.motif_12_model ...
              tblMotifs.motif_13_model tblMotifs.motif_14_model ...
              tblMotifs.motif_15_model tblMotifs.motif_16_model];
tbl.pMotifRnd = [tblMotifs.motif_1_random tblMotifs.motif_2_random ...
              tblMotifs.motif_3_random tblMotifs.motif_4_random ...
              tblMotifs.motif_5_random tblMotifs.motif_6_random ...
              tblMotifs.motif_7_random tblMotifs.motif_8_random ...
              tblMotifs.motif_9_random tblMotifs.motif_10_random ...
              tblMotifs.motif_11_random tblMotifs.motif_12_random ...
              tblMotifs.motif_13_random tblMotifs.motif_14_random ...
              tblMotifs.motif_15_random tblMotifs.motif_16_random];
clear tblMotifs; 
           
%%  
fprintf('--- samples = %d ---\n',sample);
dev = nan(2,16,length(layers)); 
for layerID = layers

    str = ['H01-' prefix '-L' num2str(layerID) '-'];
    idx = strcmp(tbl.A,[str num2str(sample)]); 
    idxModel = strcmp(tbl.A,[str num2str(sample) '#null-model']);

    if sum(idx)==0 || sum(idxModel)==0
       error('sample %s not found!',[str num2str(sample)]); 
    end
    idxtmp = (tbl.pMotif(idx,:)==0);
    fprintf('L%d: %d valid motifs\n',layerID,sum(~idxtmp));
    dev_emp = tbl.pMotif(idx,:)./tbl.pMotifRnd(idx,:);
    dev_model = tbl.pMotif(idxModel,:)./tbl.pMotifRnd(idxModel,:);
    dev_emp = dev_emp';
    dev_model = dev_model';
    %dev_model(idxtmp) = nan;
    dev_emp(idxtmp) = nan; 

    dev(1,:,layerID==layers) = dev_emp;
    dev(2,:,layerID==layers) = dev_model;
end

dev_mean = mean(dev,3,'omitnan'); 
dev_sd = std(dev,[],3,'omitnan'); 

%%
f1 = figure(1);
clf;
b = bar([1:16],dev_mean,1,'EdgeColor','none', ...
        'FaceColor','k','BaseValue',1); 
b(1).FaceColor = 'g';
b(2).FaceColor = 'k';
hold on;
offset1 = 0.15; 
offset2 = 0.13; 
lenBar = 0.12; 
for i = 1:16
    if dev_mean(1,i)>1
        plot([i i]-offset2,dev_mean(1,i)+[0 dev_sd(1,i)],'g-');
        plot(i-offset2+[-lenBar lenBar],[1 1].*dev_mean(1,i)+dev_sd(1,i),'g-');
    else
        plot([i i]-offset2,dev_mean(1,i)-[0 dev_sd(1,i)],'g-');
        plot(i-offset2+[-lenBar lenBar],[1 1].*dev_mean(1,i)-dev_sd(1,i),'g-');
    end
        
    if dev_mean(2,i)>1
        plot([i i]+offset1,dev_mean(2,i)+[0 dev_sd(2,i)],'k-'); 
        plot(i+offset1+[-lenBar lenBar],[1 1].*dev_mean(2,i)+dev_sd(2,i),'k-');
    else
        plot([i i]+offset1,dev_mean(2,i)-[0 dev_sd(2,i)],'k-'); 
        plot(i+offset1+[-lenBar lenBar],[1 1].*dev_mean(2,i)-dev_sd(2,i),'k-');
    end
end
set(gca,'TickDir','out','Box','off','XLim',[0.5 15.5],'XTick', ...
    [1:16],'YScale','log','YLim',[0.1 2e8]); 
ylabel('deviation');

set(f1,'PaperPositionMode','auto','Position',[0 0 500 350]);
print(f1,'-dsvg','-r600',[figurePath 'Motifs_' prefix '_' num2str(sample) '.svg']);
% saveas(f1,[figurePath 'Motifs_' prefix '_' num2str(sample)  '.jpg']);