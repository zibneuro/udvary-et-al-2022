%% Display results of permutation test of correlation values between 
% empirical and predicted connection probabilities
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';

%% Correlation
load([matlabPath 'data\correlationValuesRandom'], ...
            'r','rRandom','numTrails','rej');
        
rBin = [-1:0.05:1];

fprintf('---------\n');
fprintf('RandomPermutationTest: Correlation\n');
fprintf('Rejections (#trails = %d): %d (p = %.1e)\n',numTrails,rej,rej/numTrails);
fprintf(' RndCorr: M = %.2f SD = %.2f Range = [%.2f %.2f]\n', ...
            mean(rRandom),std(rRandom),min(rRandom),max(rRandom)); 
fprintf(' ModelCorr: r = %.2f\n',r); 
% Figure
f1 = figure(1);
clf;
h1 = histogram(rRandom,rBin,'FaceColor','k','EdgeColor','none','FaceAlpha',1, ...
    'Normalization','Probability');
hold on;
plot([r r],[0 max(h1.Values)],'b-');
xlabel('correlation'); 
set(gca,'YLim',[0 max(h1.Values)],'XLim',[-1 1],'Box','off','TickDir','out',...
    'YTick',[]);

% set(f1,'PaperPositionMode', 'auto','Position',[0 0 200 100]);
% print(f1,'-dsvg','-r600',[figuresPath 'RandomPermutationTest_Results.svg']);