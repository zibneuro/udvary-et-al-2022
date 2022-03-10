%% Figure 4G
% sparsity vs. heterogeneity/correlation (lambda) as colored image
% plus triplet motif deviaiton from barrel cortex model as dots
% 
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
idxExamples = [56 84 120]; % L5IT L5PT L6ACC

p_avg = tbl.avg_all.*100;
p_sd = tbl.sd_all;
p_cv = tbl.sd_all./tbl.avg_all; 
motifIDs = [13 10 1];
edgeList = [6 5 4 4 4 3 3 4 3 3 2 2 2 2 1];
col = 'rgb';

%%
lambdaMapped = nan(size(p_avg)); 
for i = 1:length(p_avg)
    [~,l] = getTheoremPrediction([matlabPath ...
        'data\math_model\Simulation_16_02_2020_19_23'], ...
        p_avg(i)/100,p_sd(i)^2); 
    lambdaMapped(i) = l; 
end

ylims = [0 40]; 
xlims = [0 1]; % lambda
if min(p_avg)<ylims(1) || max(p_avg)>ylims(2) ...
        || min(lambdaMapped)<xlims(1) || max(lambdaMapped)>xlims(2)
    error('Values outside of axis limits!');
end

%% Simulation results of Theorem 
load([matlabPath 'data\math_model\Simulation_16_02_2020_19_23'], ... 
    'lambdaList','mu_sample','devMotif','lambda'); 
numGridValues = 20; 
lambdaBins = linspace(lambdaList(1),lambdaList(end),numGridValues); % Lambda
meanBins = linspace(0,1,numGridValues); % Mean values (use as many as for lambda)
lambdaDelta = lambdaBins(2)-lambdaBins(1);
meanDelta = meanBins(2)-meanBins(1);
mu = mean(mu_sample,3);
[N,~,~,binL,binM] = histcounts2(lambda(:,1:end-1),mu(:,1:end-1), ...
                        lambdaBins,meanBins);

%%
f1 = figure(1);
clf;
c = 1; 
for motifID = motifIDs
    
    % Barrel Cortex // Deviation    
    dev_tmp = tbl.pMotif(:,motifID)./tbl.pMotifRnd(:,motifID);
%     devColor = getBlueRedColorMatrix(dev_tmp,1);
    devColor = getBlueRedColorMatrix(log10(dev_tmp),0);
    devColor = squeeze(devColor); 
    [~,idxSort] = sort(dev_tmp); 
    fprintf('%d: Range of deviation = [%.2f %.2f]\n', ...
                    motifID,min(dev_tmp(:)),max(dev_tmp(:))); 
          
    % Theory // Edges
    k = edgeList(motifID); 
    % Note skip lambda = 1 (because pMotif = 0 sometimes -> dev = 0); 
    dev_theory = squeeze(devMotif(:,1:end-1,k+1));
    N2 = zeros(size(N)); 
    C2 = zeros(size(N));

    for x = 1:numel(binL)            
        N2(binL(x),binM(x)) = N2(binL(x),binM(x)) + dev_theory(x);
        C2(binL(x),binM(x)) = C2(binL(x),binM(x)) + 1; 
    end
    Nnorm = N2./C2; 
    
    % Crop only to 40%
    %Nnorm(:,meanBins>ylims(2)/100) = nan;   
    devColor_theory = getBlueRedColorMatrix(log10(Nnorm'),0);
    
%     fprintf([' >> %d Edges: Range of deviation in each grid point ' ...
%                 '= [%.2e %.2e]\n'],k,min(dev_theory(:)),max(dev_theory(:))); 
    fprintf([' >> %d Edges: Range of deviation in each grid point ' ...
                '= [%.2e %.2e]\n'],k,nanmin(Nnorm(:)),nanmax(Nnorm(:)));
            
    % Plot
    subplot(3,1,c);
    hold on;
    imagesc(lambdaBins(1:end-1)+lambdaDelta/2, ...
                (meanBins(1:end-1)+meanDelta/2).*100,devColor_theory); 
    
    scatter(lambdaMapped(idxSort),p_avg(idxSort),20, ...
                    devColor(idxSort,:),'o','filled','MarkerEdgeColor','none');
    for k = 1:length(idxExamples)
        scatter(lambdaMapped(idxExamples(k)),p_avg(idxExamples(k)),20, ...
                    devColor(idxExamples(k),:),'filled', ...
                    'MarkerEdgeColor',col(k));
    end
    
    ylim(ylims);
    xlim(xlims);
    set(gca,'TickDir','out','Box','off','YTick',[0:20:40],'XTick',[0:0.5:2]);
%     xlabel('lambda(p)');
%     ylabel('mean(p)'); 
    title([num2str(motifID) ' (#edges=' num2str(k) ')']);
    axis square;
    c = c +1; 
    
%     [r,p] = corr((p_avg)',dev_tmp'); 
%     fprintf('CORR(AVG,DeviationBCModel): r = %.2f; p = %.2e\n',r,p);
%     [r,p] = corr((p_cv)',dev_tmp'); 
%     fprintf('CORR(CV,DeviationBCModel): r = %.2f; p = %.2e\n',r,p);
    tblResult = table(tbl.A,tbl.B,tbl.C,p_avg,lambdaMapped,dev_tmp, ...
            'VariableNames',{'Cell_A','Cell_B','Cell_C', ...
            'Mean_ConnectionProbability','Lambda','deviation'});
    writetable(tblResult,[tablePath 'MathModel_vs_TripletMotif_' ...
            num2str(motifID) '.csv']);
end

%% Save everything as csv file
set(f1,'PaperPositionMode', 'auto','Position',[0 0 500 900]);
print(f1,'-dsvg','-r600',[figurePath 'DevTripletCombinations_All_v2.svg']);