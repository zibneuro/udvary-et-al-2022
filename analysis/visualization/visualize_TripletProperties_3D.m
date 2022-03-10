%% Fig 4F
% 3D: mean vs CV vs correlation vs dev (colorcoded)
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

writetable(tbl,[tablePath 'CorrelationInDegreesForTriplets.csv']);
motifIDs = [13 10 1];

%%
p_avg = tbl.avg_all.*100;
p_cv = tbl.sd_all./tbl.avg_all; 
Rmapped = tbl.R_avg;
fprintf('P_AVG = [%.2f %.2f]\n',min(p_avg),max(p_avg));
fprintf('P_CV = [%.2f %.2f]\n',min(p_cv),max(p_cv));
fprintf('R = [%.2f %.2f]\n',nanmin(Rmapped),nanmax(Rmapped));
fprintf('n = %d\n',numel(p_cv));

xlims = [0 40];
ylims = [0.4 2.2];
zlims = [-0.25 0.9];

if min(p_avg)<xlims(1) || max(p_avg)>xlims(2) ...
        || min(p_cv)<ylims(1) || max(p_cv)>ylims(2) ...
        || nanmin(Rmapped)<zlims(1) || nanmax(Rmapped)>zlims(2)
    error('Values outside of axis limits!');
end

%%
f1 = figure(1);
clf;
c = 1; 
for motifID = motifIDs
    
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
       error('something is wrong here'); 
    end
    
%     devColor = getBlueRedColorMatrix(dev_tmp,1);
    devColor = getBlueRedColorMatrix(log10(dev_tmp),0);
    devColor = squeeze(devColor); 
    [~,idxSort] = sort(dev_tmp); 
    fprintf('%d: Range of deviation = [%.2f %.2f]\n', ...
                    motifID,min(dev_tmp(:)),max(dev_tmp(:))); 

    figure(1);
    subplot(1,3,c); 
    scatter3(p_avg(idxSort),p_cv(idxSort),Rmapped(idxSort),20, ...
                    devColor(idxSort,:),'filled','MarkerEdgeColor','none');
    hold on;
%     for c = 1:length(idxExamples)
%         scatter3(p_avg(idxExamples(c)),p_cv(idxExamples(c)), ...
%                     Rmapped(idxExamples(c)),20, ...
%                     devColor(idxExamples(c),:),'filled', ...
%                     'MarkerEdgeColor',col(c));
%     end
    
    set(gca,'TickDir','out','Box','off', ...
            'ZLim',zlims,'YLim',ylims,'XLim',xlims, ...
            'ZTick',[-0.3:0.3:0.9],'YTick',[0:1:4],'XTick',[0:20:40]);
    xlabel('M');
    ylabel('CV'); 
    zlabel('R');
    title(['Motif ' num2str(motifID)]);
    axis square;
    c = c+1; 
    
    tblResult = table(tbl.A,tbl.B,tbl.C,p_avg,p_cv,Rmapped,pMotif_tmp, ...
            pMotifRnd_tmp,dev_tmp, ...
            'VariableNames',{'Cell_A','Cell_B','Cell_C', ...
            'Mean_ConnectionProbability','CV_ConnectionProbability', ...
            'Correlation_R','p(model)','p(random)','deviation'});
    writetable(tblResult,[tablePath 'TripletProperties_3D_Motif_' ...
            num2str(motifID) '.csv']);
end

%%
set(f1,'PaperPositionMode', 'auto','Position',[0 0 900 500]);
print(f1,'-painters','-dsvg','-r600',[figurePath ...
            'TripletProperties_3D.svg']);