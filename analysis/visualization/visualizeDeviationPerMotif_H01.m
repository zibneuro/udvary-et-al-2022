%% Figure 7D
% predicted vs. observed motif deviations per Motif and subsample
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc
warning on;

matlabPath = 'D:\udvary-et-al-2022\analysis\';
figurePath = [matlabPath 'output\Figures\']; 
prefix = 'PYR';
datapath = [matlabPath 'data\H01\motifs_' prefix '\h01-subpopulations\'];
tblMotifs = readtable([datapath 'probabilities_16_h01-subpopulations.csv']);
col = [107 62 152; 246 153 153; 27 175 72; 250 162 45; 31 120 180]./255;
layers = 2:6; 
samples = [50 100 150 200]; 
symbs = 'osd^';

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
f1 = figure(1);
clf; 

targetMotif = 11; 
pMotif = tbl.pMotif(:,targetMotif);
pMotifRnd = tbl.pMotifRnd(:,targetMotif);
dev = pMotif./pMotifRnd; 
label = cell(0); 
devtmp = [];

for layerID = layers
    str = ['H01-' prefix '-L' num2str(layerID) '-'];

    for s = samples
        label{end+1} = [str num2str(s)];  
        idx = strcmp(tbl.A,[str num2str(s)]); 
        idxModel = strcmp(tbl.A,[str num2str(s) '#null-model']);

        if sum(idx)==0 || sum(idxModel)==0
           error('sample %s not found!',[str num2str(s)]); 
        end

        devtmp = [devtmp; dev(idx) dev(idxModel)]; 
    end
end
devtmp = devtmp';
idxtmp = devtmp(1,:)==0;
devtmp(:,idxtmp) = nan; 
devLog = log10(devtmp);
[r,p] = corr(devLog(1,~idxtmp)',devLog(2,~idxtmp)');
fprintf('Motif %d: R = %.2f p = %.2e\n',targetMotif,r,p);

figure(1);
subplot(1,2,1); 
plot([min(devLog(:)) max(devLog(:))],[min(devLog(:)) max(devLog(:))], ...
        '--','Color',[0.5 0.5 0.5]);
hold on; 
for layerID = layers
    for s = samples
        idxtmp = contains(label,['L' num2str(layerID) '-' num2str(s)]); 
        plot(devLog(2,idxtmp),devLog(1,idxtmp),symbs(s==samples), ...
            'MarkerEdgeColor','none','MarkerFaceColor',col(layerID==layers,:));
    end
end
set(gca,'TickDir','out','Box','off','XLim',[-0.5 0.5],'YLim',[-1 0.5]);
xlabel('log(dev) prediction');
ylabel('log(dev) empirical');
axis square; 

%%
devtmp = [];
label = cell(0); 

for motifIDX = [1:15]
    
    if motifIDX==targetMotif
        continue;
    end
    
    pMotif = tbl.pMotif(:,motifIDX);
    pMotifRnd = tbl.pMotifRnd(:,motifIDX);
    dev = pMotif./pMotifRnd; 
    
    for layerID = layers
       str = ['H01-' prefix '-L' num2str(layerID) '-'];
                
        for s = samples
            label{end+1} = [str num2str(s)];  
            idx = strcmp(tbl.A,[str num2str(s)]); 
            idxModel = strcmp(tbl.A,[str num2str(s) '#null-model']);
            
            if sum(idx)==0 || sum(idxModel)==0
               error('sample %s not found!',[str num2str(s)]); 
            end
            devtmp = [devtmp; dev(idx) dev(idxModel)]; 
        end
    end
end

devtmp = devtmp';
idxtmp = devtmp(1,:)==0;
devtmp(:,idxtmp) = nan; 
devLog = log10(devtmp);
[r,p] = corr(devLog(1,~idxtmp)',devLog(2,~idxtmp)');
fprintf('Other: R = %.2f p = %.2e\n',r,p);

subplot(1,2,2); 
plot([min(devLog(:)) max(devLog(:))],[min(devLog(:)) max(devLog(:))], ...
        '--','Color',[0.5 0.5 0.5]);
hold on;
for layerID = layers
    for s = samples
        idxtmp = contains(label,['L' num2str(layerID) '-' num2str(s)]); 
        plot(devLog(2,idxtmp),devLog(1,idxtmp),symbs(s==samples), ...
            'MarkerEdgeColor','none','MarkerFaceColor',col(layerID==layers,:));
    end
end

set(gca,'TickDir','out','Box','off'); 
xlabel('log(dev) prediction');
ylabel('log(dev) empirical');
axis square;

set(f1,'PaperPositionMode','auto','Position',[0 0 500 300]);
% saveas(f1,[figurePath 'Motif' num2str(targetMotif) '_vs_Rest_' prefix '.jpg']);
print(f1,'-dsvg','-r600',[figurePath 'Motif' num2str(targetMotif) ...
            '_vs_Rest_' prefix '.svg']);