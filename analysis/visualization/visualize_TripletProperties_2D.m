%% Figure 3F
% sparsity vs. heterogeneity (CV)
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
figurePath = [matlabPath 'output\Figures\']; 
tablePath = [matlabPath 'output\Tables\']; 
load([matlabPath 'data\TripletMotifs_C2.mat'],'tbl'); 

% Delete all INH and VPM
idxDel = strcmp(tbl.A,'VPM') | strcmp(tbl.B,'VPM') | strcmp(tbl.C,'VPM') ...
    | strcmp(tbl.A,'INH') | strcmp(tbl.B,'INH') | strcmp(tbl.C,'INH');
tbl(idxDel,:) = []; 

idxExamples = [1 84 171]; % 1: L2PY-L2PY-L2PY; 84: L5PT-L5PT-L5PT
% 171: L6CT-L4PY-L4PY
col = 'rgb';

%%
p_avg = tbl.avg_all.*100;
p_cv = tbl.sd_all./tbl.avg_all; 
fprintf('P_AVG = [%.2f %.2f]\n',min(p_avg),max(p_avg));
fprintf('P_CV = [%.2f %.2f]\n',min(p_cv),max(p_cv));
fprintf('n = %d\n',numel(p_cv));

xlims = [0 40];
ylims = [0 2.5];

if min(p_avg)<xlims(1) || max(p_avg)>xlims(2) ...
        || min(p_cv)<ylims(1) || max(p_cv)>ylims(2)
    error('Values outside of axis limits!');
end

m1 = max(abs(sum(tbl.pMotif,2)-1));
m2 = max(abs(sum(tbl.pMotifRnd,2)-1));
fprintf('Max deviation from 1: %.2e %.2e\n',m1,m2);

motifID = 10;
pMotif_tmp = tbl.pMotif(:,motifID);
pMotifRnd_tmp = tbl.pMotifRnd(:,motifID); 
dev_tmp = pMotif_tmp./pMotifRnd_tmp;

% Check cases where pMotifRnd and pMotif is 0
idx01 = (pMotifRnd_tmp==0);
idx02 = (pMotif_tmp==0);
% if both zero, no deviation
dev_tmp(idx01 & idx02) = 1; 
idxtmp = xor(idx01,idx02); 
if sum(idxtmp)>0
   warning('something is wrong here'); 
end

%     devColor = getBlueRedColorMatrix(dev_tmp,1);
devColor = getBlueRedColorMatrix(log10(dev_tmp),0);
devColor = squeeze(devColor); 
[~,idxSort] = sort(dev_tmp,'descend'); 
fprintf('%d: Range of deviation = [%.2f %.2f]\n', ...
                motifID,min(dev_tmp(:)),max(dev_tmp(:))); 

[r1,p1] = corr(p_avg,dev_tmp); 
fprintf('CORR(AVG,DeviationBCModel): r = %.2f; p = %.2e\n',r1,p1);
[r2,p2] = corr(p_cv(p_avg>0),dev_tmp(p_avg>0)); 
fprintf('CORR(CV,DeviationBCModel): r = %.2f; p = %.2e\n',r2,p2);
            
%%
f1 = figure(1);
clf; 
scatter(p_avg(idxSort),p_cv(idxSort),20, ...
                devColor(idxSort,:),'filled','MarkerEdgeColor','none'); 
hold on;
for c = 1:length(idxExamples)
    scatter(p_avg(idxExamples(c)),p_cv(idxExamples(c)),20, ...
                devColor(idxExamples(c),:),'filled','MarkerEdgeColor',col(c));
end
set(gca,'TickDir','out','Box','off','YLim',ylims,'XLim',xlims, ...
        'YTick',[0:1:8],'XTick',[0:20:100]);
xlabel('M');
ylabel('CV'); 
title(['Motif ' num2str(motifID)]);
axis square;
set(f1,'PaperPositionMode', 'auto','Position',[0 0 200 200]);
print(f1,'-painters','-dsvg','-r600',[figurePath ...
            'TripletProperties_2D.svg']);
        
%%
tblResult = table(tbl.A,tbl.B,tbl.C,p_avg,p_cv,pMotif_tmp,pMotifRnd_tmp,dev_tmp, ...
            'VariableNames',{'Cell_A','Cell_B','Cell_C', ...
            'Mean_ConnectionProbability','CV_ConnectionProbability', ...
            'p(model)','p(random)','deviation'});
writetable(tblResult,[tablePath 'TripletProperties_2D_Example.csv']);