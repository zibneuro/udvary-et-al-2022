%% Figure S5C
% Red-Blue Image Density plot of motif deviation accroding to mathematical
% model with correlated connectivity
%
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
addpath(genpath([matlabPath 'functions\']));
addpath(genpath([matlabPath 'math_model_corr_connectivity\']));
figurePath = [matlabPath 'output\Figures\']; 

%% Figure
load([matlabPath 'data\math_model\Simulation_16_02_2020_19_23']); 
numGridValues = 20; 
lambdaBins = linspace(lambdaList(1),lambdaList(end),numGridValues); % Lambda
meanBins = linspace(0,1,numGridValues); % Mean values (use as many as for lambda)
lambdaDelta = lambdaBins(2)-lambdaBins(1);
meanDelta = meanBins(2)-meanBins(1);
mu = mean(mu_sample,3);
[N,~,~,binL,binM] = histcounts2(lambda(:,1:end-1),mu(:,1:end-1), ...
                        lambdaBins,meanBins);
edgeList = [1:6]; 

f1 = figure(1);
clf;
c = 1; 
for k = edgeList

    % Note skip lambda = 1 (because pMotif = 0 sometimes -> dev = 0); 
    dev = squeeze(devMotif(:,1:end-1,k+1));
    N2 = zeros(size(N)); 
    C2 = zeros(size(N));

    for x = 1:numel(binL)            
        N2(binL(x),binM(x)) = N2(binL(x),binM(x)) + dev(x);
        C2(binL(x),binM(x)) = C2(binL(x),binM(x)) + 1; 
    end

    fprintf([' %d Edges: Range of deviation in each grid point ' ...
                '= [%.2e %.2e]\n'],k,min(dev(:)),max(dev(:))); 

    Nnorm = N2./C2; 

    % Plot in grid
    devColor = getBlueRedColorMatrix(log10(Nnorm'),0);

    subplot(2,3,c);
    c = c + 1; 
    imagesc(lambdaBins(1:end-1)+lambdaDelta/2, ...
                meanBins(1:end-1)+meanDelta/2,devColor);     
    set(gca,'TickDir','out','Box','off','YTick',[0:0.2:1], ...
            'XTick',[0:0.2:1],'YDir','normal','Xlim',[0 1],'YLim',[0 1]);
    %xlabel('\lambda');
    ylabel('mean(p)'); 
    title([num2str(k) ' Edges; dev = [' num2str(min(dev(:)),'%.2f') ' ' ...
                num2str(max(dev(:)),'%.2e') ']']);
    axis square;
end

set(f1,'PaperPositionMode', 'auto','Position',[0 0 700 350]);
print(f1,'-dsvg','-r600',[figurePath 'TheoremPredictionLambda.svg']);