%% Figure 4D (top)
% Histogram of L6CC to L5PT; L5IT to L5PT; L5PT to L5PT
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
figurePath = [matlabPath 'output\Figures\']; 
tablePath = [matlabPath 'output\Tables\']; 

binsz = 0.02; 
XSample = [0:binsz:1+binsz];

%%
f1 = figure(1);
clf;

maxVal = 0; 

for ii = 1:3
    
    switch ii
        case 1
            preType = 'L5IT';
            postType = 'L5PT';
            col = [0 1 0];
        case 2
            preType = 'L5PT';
            postType = 'L5PT';  
            col = [1 0.5 0];
        case 3
            preType = 'L6ACC';
            postType = 'L5PT';  
            col = [0 0 1];
    end

    load([matlabPath 'data\CellMatrix_C2\' preType '_' postType ...
                        '.mat']);
    
    % Load Innervation values
    p = 1-exp(-I.I(:));
    N = histcounts(p,XSample);
    maxVal = max([maxVal N]); 
    
    fprintf('%s-%s: M = %.2f SD = %.2f CV = %.2f\n',preType,postType, ...
            mean(p),std(p),std(p)/mean(p));
    
    % Plot
    hold on;
    stairs([XSample(1) XSample(1:end-1) XSample(end)].*100,[0 N 0], ...
                'Color',col); 
    ylabel('occurrences');
    xlabel(['P(A,B) (%)']);
    set(gca,'TickDir','out','Box','off','XLim',[-2 100],'XTick',[0:25:100], ...
        'YLim',[0 maxVal],'YTick',[0 maxVal/2 maxVal],'YTickLabel',''); 
    
    tbl = table();
    tbl.connection_probability = round(p,4); 
    writetable(tbl,[tablePath preType '_' postType '_connection_probabilities.csv']);
end

%% Save figure
set(f1,'PaperPositionMode','auto','Position',[0 0 250 150]); 
print(f1,'-dsvg','-r600',[figurePath 'ExampleDistributions_CellType_EmpDist.svg']); 