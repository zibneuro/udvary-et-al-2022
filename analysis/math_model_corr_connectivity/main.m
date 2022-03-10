%% Numerical simulation of mathematical model of correlated connectivity
% Script computes the deviation of observing seven triplet motifs 
% (connected by 0, 1, 2, 3, 4, 5 or 6 edges) in a correlated network model 
% vs. observing the seven triplet motifs in a random network 
% (only constraint by the mean connection probability). 
% Deviation with respect to motif probabilities. 
%
% Script saves results of simulation in local directory
% Script generates two figures of results
% Figure 1: Deviation of each of the seven motifs
% Figure 2: Only significant deviation of each of the seven motifs
% 
% Model developed by Jakob H. Macke (University Tuebingen)
% Model implemented by Daniel Udvary (Max Planck Institute for 
% Neurobiology of Behavior – caesar)
% Date: Feb, 9th, 2022
% 
% Relevant literature: 
% - Model reported in:
%   The Impact of Neuronal Structure on Cortical Network Architecture
%   Daniel Udvary, Philipp Harth, Jakob H. Macke, Hans-Christian Hege, 
%   Christiaan P.J. de Kock, Bert Sakmann, Marcel Oberlaender
%   Cell Reports (2022)
%
% - Model closely related to a model studied in: 
%   J. H. Macke, M. Opper, M. Bethge, 
%   Common input explains higher-order correlations and 
%   entropy in a simple model of neural population activity. 
%   Physical Review Letters 106, 208102 (2011).

clear all
close all
clc

%%
numSamples = 1e5;
numTrials = 10;
numEdges = 6; 

% Number of simulated gamma values (gamma: degree of connectivity)
numGammaValues = 250; 
gammaValues = linspace(-2,2,numGammaValues);

% Number of simulated lambda values (lambda: measure of correlation 
% and variance; degree of correlation between S (shared source) and 
% Ti (private source)
numLambdaValues = 250; 
lambdaValues = linspace(0,1,numLambdaValues); 

% Get norm cdf (note: Input is standard deviation)
L = @(x,mu,sd) 1 - normcdf(x,mu,sd);

% Initalize variables
mu_analytical = nan(length(gammaValues),length(lambdaValues)); 
mu_sample = nan(length(gammaValues),length(lambdaValues),numTrials); 
var_sample = nan(size(mu_sample)); 
probMotif = nan(length(gammaValues),length(lambdaValues),numEdges+1); 
probRandom = nan(length(gammaValues),length(lambdaValues),numEdges+1); 
devMotif = nan(length(gammaValues),length(lambdaValues),numEdges+1); 
minZValues = nan(length(gammaValues),length(lambdaValues),numEdges+1);  
lambda = nan(size(mu_analytical)); 
gamma = nan(size(mu_analytical)); 

% Threshold for zscore (to determine significant deviations)
zThreshold = 4; 

%%
fprintf('Start Simulation\n-------------\n');
rng(0);
simulateResults = 0; 
sum_imprecision = 0;
mean_imprecision = 0;

for i = 1:length(gammaValues)
    
    fprintf('%d / %d (gamma = %.2f)\n',i,length(gammaValues),gammaValues(i));
    
    % Degree of connectivity
    gamma_tmp = gammaValues(i);
    
    for j = 1:length(lambdaValues)
        
        gamma(i,j) = gamma_tmp; 
        
        % Variance of shared source S
        % -> degree of correlation
        % -> large lambda, large contribution of shared source to Z
        % -> higher correlation
        lambda_var = lambdaValues(j); 
        lambda(i,j) = lambda_var;
        
        % Variance of private source Ti
        n_var = 1 - lambda_var; 
        
        if n_var<0
            error('Something is wrong! Variance is negative!');
        end
        
        mu_analytical(i,j) = L(0,gamma_tmp,sqrt(lambda_var+n_var));        
        mu_tmp = nan(1,numTrials);
        var_tmp = nan(1,numTrials);
        pMotif_tmp = nan(numEdges+1,numTrials);
        pRandom_tmp = nan(numEdges+1,numTrials);
        
        for t = 1:numTrials
            result = simulateTheorem(numSamples,gamma_tmp, ...
                                            lambda_var,numEdges, ...
                                            mu_analytical(i,j));
            
            % Update Imprecisions
            mean_imprecision = max([result.mean_imprecision ...
                                        mean_imprecision]);
            sum_imprecision = max([result.sum_imprecision ...
                                        sum_imprecision]);
                  
            % Mean and variance over connection probabilities
            mu_tmp(t) = result.p_mu;
            var_tmp(t) = result.p_var;
            
            % Motif probabilities
            pMotif_tmp(:,t) = [result.pMotif0 result.pMotif];
            pRandom_tmp(:,t) = [result.pRandom0 result.pRandom];
        end
        
        % Check whether differences between Motif probabilities 
        % are significant
        pMotif = mean(pMotif_tmp,2);
        pRandom = mean(pRandom_tmp,2); 
        se = std(pMotif_tmp,0,2)./size(pMotif_tmp,2);
        z = (pMotif - pRandom)./se;
        z(pMotif==pRandom) = 0; 
        
        devMotif(i,j,:) = pMotif./pRandom;
        minZValues(i,j,:) = z; 
        mu_sample(i,j,:) = mu_tmp; 
        var_sample(i,j,:) = var_tmp; 
        probMotif(i,j,:) = pMotif; 
        probRandom(i,j,:) = pRandom; 
    end
end

%% Store results to hard drive
save(['Simulation_' datestr(now, 'dd_mm_yyyy_HH_MM') '.mat']); 
fprintf('Simulation completed and stored\n-------------\n');

%% Display results
fprintf('Imprecision = %.2e\n',sum_imprecision);
fprintf('Range of motif deviation = [%.2f %.2f]\n', ...
            min(devMotif(:)),max(devMotif(:))); 
devMean = mu_sample./mu_analytical;
fprintf('Range of mean deviation = [%.2f %.2f]\n', ...
            min(devMean(:)),max(devMean(:))); 
fprintf('Range of z values = [%.2f %.2f]\n', ...
            min(minZValues(:)),max(minZValues(:))); 
r = sum(minZValues(:)<zThreshold)./numel(minZValues)*100;
fprintf('z Values below %.2f: %.2f%%\n',zThreshold,r);   
fprintf('-------------\n');

%% Visualize results
% Map results on grid
numGridValues = 20; 
% Lambda
lambdaBins = linspace(0,1,numGridValues);
% Mean values (use as many as for lambda)
meanBins = linspace(0,1,numGridValues); 
lambdaDelta = lambdaBins(2)-lambdaBins(1);
meanDelta = meanBins(2)-meanBins(1);

% Opt==1 Show all deviaton values 
% Opt==2 Set deviaton values with z-value<zThreshold to 1 (i.e., no deviation)
optStr = {'all',['z<' num2str(zThreshold,'%.2f')]};
for opt = 1:2
    
    f1 = figure(opt);
    clf;    
    mu = mean(mu_sample,3);
    [N,~,~,binL,binM] = histcounts2(lambda,mu,lambdaBins,meanBins);
    
    fprintf('---- %s ----\n',optStr{opt});
    
    for k = 1:numEdges+1

        dev = squeeze(devMotif(:,:,k));
        N2 = zeros(size(N)); 
        C2 = zeros(size(N));
        
        if opt==2
           z = squeeze(abs(minZValues(:,:,k))); 
           dev(z<zThreshold) = 1; 
        end
        
        for x = 1:numel(binL)            
            N2(binL(x),binM(x)) = N2(binL(x),binM(x)) + dev(x);
            C2(binL(x),binM(x)) = C2(binL(x),binM(x)) + 1; 
        end

        fprintf([' %d Edges: Range of deviation in each grid point ' ...
                    '= [%.2f %.2f]\n'],k-1,min(dev(:)),max(dev(:))); 
        fprintf('   Range of samples in each grid voxel = [%.2f %.2f]\n', ...
            min(C2(:)),max(C2(:))); 

        Nnorm = N2./C2; 

        % Plot in grid
        devColor = getBlueRedColorMatrix(log10(Nnorm'),0);
        subplot(2,4,k);
        imagesc(lambdaBins(1:end-1)+lambdaDelta/2, ...
                    meanBins(1:end-1)+meanDelta/2,devColor); 
        hold on;
        set(gca,'TickDir','out','Box','off','YTick',[0:0.2:1], ...
                'XTick',[0:0.2:1],'YDir','normal','Xlim',[0 1],'YLim',[0 1]);
        xlabel('\lambda');
        ylabel('mean(p)'); 
        title([num2str(k-1) ' Edges; dev = [' num2str(min(dev(:)),'%.2f') ...
                    ' ' num2str(max(dev(:)),'%.2e') ']']);
        axis square;
    end
end