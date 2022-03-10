%% Figure S6B
% - Boxplot of Correlation values per Cell Type Combination
%
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
figurePath = [matlabPath 'output\Figures\']; 
tablePath = [matlabPath 'output\Tables\']; 
load([matlabPath 'data\CorrData.mat']);

%% Correlation Differences for each postsynaptic cell type
CellTypeList = unique(PostType);
N = 66; 
numCT = length(CellTypeList); 
rAll = nan(N,numCT); 
%pAll = nan(N,numCT); 

for ct = 1:numCT
    idxtmp = strcmp(CellTypeList{ct},PostType);
    rAll(:,ct) = R(idxtmp);
    %pAll(:,ct) = P(idxtmp);
end

%%
f1 = figure(1); 
clf;
wdth = 0.7;
lnwdth = 2; 
plot([0 0],[0 numCT+1],'k--');
hold on; 
boxplot(rAll,'Color','k','Symbol','+','Orientation','Horizontal');

set(gca,'TickDir','out','Box','off','YLim',[0 numCT]+0.5,'YTick', ...
    [1:numCT],'YTickLabel',CellTypeList,'XLim',[-1 1],'YDir','reverse');
xlabel('correlated input onto postsynaptic cell type');
ylabel('correlation r'); 

set(f1,'PaperPositionMode', 'auto','Position',[0 0 300 500]);
print(f1,'-dsvg','-r600',[figurePath 'CorrelationInputOverview.svg']);

%% Display results
for idxCTPost = 1:length(CellTypeList)
   rtmp = rAll(:,idxCTPost);  
   fprintf('%s: M = %.2f SD = %.2f Range = [%.2f %.2f]\n', ...
            CellTypeList{idxCTPost},mean(rtmp),std(rtmp),min(rtmp),max(rtmp)); 
end

rtmp = rAll(:);  
fprintf('All: M = %.2f SD = %.2f Range = [%.2f %.2f]\n', ...
        mean(rtmp),std(rtmp),min(rtmp),max(rtmp)); 
    
%% Display min correlations
[Rsorted,sortIDX] = sort(R); 
for i = 1:10
    fprintf('%.2f (%s,%s)->%s\n',Rsorted(i), ...
        PreType1{sortIDX(i)},PreType2{sortIDX(i)},PostType{sortIDX(i)});
end

%% Results as csv
tbl = table();
tbl.PreType1 = PreType1';
tbl.PreType2 = PreType2'; 
tbl.PostType = PostType'; 
tbl.R = R'; 
writetable(tbl,[tablePath 'InDegree_Correlations_CellTypeCombinations.csv']);