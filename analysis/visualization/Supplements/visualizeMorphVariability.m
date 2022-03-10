%% Figure S2G
% sample size of morphologies vs. change in path length and contributing
% cells
%
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
figuresPath = [matlabPath 'output\Figures\']; 
tablesPath = [matlabPath 'output\Tables\']; 
numVoxels = 512; 
load([matlabPath 'data\variability\MergedAxonDend.mat'], ...
        'CV_CC','CV_Len'); 
d = load([matlabPath 'data\variability\MergedAxonDend_DiffSoma.mat'], ...
        'CV_CC','CV_Len'); 
scale = 'log';

%%
f1 = figure(1);
clf;
strLabel = {'Length','#Cells'};
tbl = table(); 
tbl.SampleSize = [1:13]';

for k = 1:length(strLabel)

    switch k
        case 1
            y = CV_Len;
            y_soma = d.CV_Len;
            tbl.CV_Morph_Length = CV_Len';
            tbl.CV_Cellular_Length = d.CV_Len';
        case 2
            y = CV_CC;
            y_soma = d.CV_CC; 
            tbl.CV_Morph_ContributingCells = CV_CC';
            tbl.CV_Cellular_ContributingCells = d.CV_CC';
    end
    
    minVal = min(y)/2; 
    x = 1:length(y); 
    
    subplot(1,2,k);
    hold on;
    plot(x,y_soma,'-o','MarkerEdgeColor','none', ...
            'MarkerFaceColor',[0.5 0.5 0.5],'Color',[0.5 0.5 0.5]);
    plot(x,y,'-o','MarkerEdgeColor','none', ...
            'MarkerFaceColor','k','Color','k');
    set(gca,'Box','off','TickDir','out','XLim',[0 length(y)]+0.5, ...
            'XTick',x,'Ylim',[minVal 1],'Yscale',scale);
    xlabel('Samples');
    ylabel(['CV ' strLabel{k}]);   
end

%% Save figure
set(f1,'PaperPositionMode','auto','Position',[0 0 600 150]); 
print(f1,'-painters','-dsvg','-r600',[figuresPath ...
        'Variability_MorphSample.svg']); 
writetable(tbl,[tablesPath 'Robustness\Variability_MorphSampleUpscaling.csv']);    