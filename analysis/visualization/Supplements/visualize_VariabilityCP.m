%% Figure S6A
% mean of CP vs. CV of CP
%
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

%%
matlabPath = 'D:\udvary-et-al-2022\analysis\';
addpath(genpath([matlabPath 'functions\']));
dataPath = [matlabPath 'data/dataParam/'];
figurePath = [matlabPath 'output\Figures\']; 
tablePath = [matlabPath 'output\Tables\']; 

%%
PreCellTypes = generateCellTypeList('v2');
PostCellTypes = PreCellTypes(1:end-1);

m = nan(length(PreCellTypes),length(PostCellTypes)); 
sd = nan(length(PreCellTypes),length(PostCellTypes)); 
skw = nan(length(PreCellTypes),length(PostCellTypes)); 

fid = fopen([tablePath 'SummaryStats_CellType_CP.csv'], 'w+');
if fid==-1
   error('Cannot write .csv file'); 
end
fprintf(fid,'PreCellType,PostCellType,');
fprintf(fid,['Mean_ConnectionProbability,SD_ConnectionProbability,' ...
            'CV_ConnectionProbability,\n']);

for preCT = 1:length(PreCellTypes)
    for postCT = 1:length(PostCellTypes)
        load([dataPath PreCellTypes{preCT} '_' ...
                        PostCellTypes{postCT} '.mat']);
                    
        m(preCT,postCT) = param_p.m;
        sd(preCT,postCT) = param_p.SD;
        skw(preCT,postCT) = param_p.skew;
        
        if param_p.skew<0
            fprintf('%.2f: %s %s\n',param_p.skew, ...
                        PreCellTypes{preCT},PostCellTypes{postCT});
        end
        
        fprintf(fid,'%s,',PreCellTypes{preCT},PostCellTypes{postCT});
        fprintf(fid,'%.6f,',m(preCT,postCT).*100,sd(preCT,postCT).*100, ...
                        sd(preCT,postCT)./m(preCT,postCT));
        fprintf(fid,'\n'); 
    end
end
fclose(fid);
cv = sd./m; 

%% scatter plot: color-code
f1 = figure(1);
clf;
scatter(cv(:),m(:),30,[0 0 0],'filled','MarkerEdgeColor','none');
set(gca,'Box','off','TickDir','out','XTick',[0:0.5:10],'YTick',[0:0.1:0.5], ...
    'XLim',[0 10],'YLim',[0 0.4]);
xlabel('CV');
ylabel('M');

set(f1,'PaperPositionMode', 'auto','Position',[0 0 800 250]);
print(f1,'-dsvg','-r600',[figurePath 'Variability_CP.svg']);

%%
fprintf('n = %d\n',numel(skw));
fprintf('Skewness = [%.2f %.2f]\n',min(skw(:)),max(skw(:)));
fprintf('Skewness<=0: %d\n',sum(skw(:)<=0));
fprintf('M>SD: %d\n',sum(m(:)>sd(:)));
fprintf('M<=SD: %d\n',sum(m(:)<=sd(:))); 
fprintf('--------------\n');
fprintf('M = %.4f %.4f [%.4f %.4f]\n',mean(m(:)),std(m(:)),min(m(:)),max(m(:)));
fprintf('SD = %.4f %.4f [%.4f %.4f]\n',mean(sd(:)),std(sd(:)),min(sd(:)),max(sd(:)));
fprintf('CV = %.4f %.4f [%.4f %.4f]\n',mean(cv(:)),std(cv(:)),min(cv(:)),max(cv(:)));
fprintf('Skew = %.4f %.4f [%.4f %.4f]\n',mean(skw(:)),std(skw(:)),min(skw(:)),max(skw(:)));