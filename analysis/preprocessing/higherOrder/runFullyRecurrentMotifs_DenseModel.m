% Run motif analysis for fully recurrent motifs of 3 to 10 nodes
% Requires: C2.mat (generated from extractC2matrix.m)
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'Z:\Documents\papers\Axon2\NullHypothesis\Final_V2\ScriptsFinal\';
inputPath = [matlabPath 'preprocessing/data/']; 
outputPath = [matlabPath 'data/higherOrder/']; 
load([inputPath 'C2.mat'],'I'); 
postfix = 'C2';
pTotal = 1-exp(-I);
clear I; 
numTrials = 1e7;
numNodesListFullyRec = [3:1:10]; 
n = size(pTotal,1); 

%%
rng(154581); 
pRec = nan(size(numNodesListFullyRec)); 
pRecRandom = nan(size(numNodesListFullyRec)); 
devFullyRec = nan(size(numNodesListFullyRec));
p_avg = nan(size(numNodesListFullyRec)); 
p_sd = nan(size(numNodesListFullyRec));

for numNodes = numNodesListFullyRec

    idDiagonal = 1:numNodes+1:(numNodes^2); % id of diagonal values
    idNonDiagonal = setdiff(1:numNodes^2,idDiagonal); 
    pRec_tmp = nan(1,numTrials); 
    
    pSelect = nan(numTrials,numNodes,numNodes); 
    
    for j = 1:numTrials
        idx = randperm(n,numNodes);
        p = pTotal(idx,idx);
        pRec_tmp(j) = prod(p(idNonDiagonal));  
        p(idDiagonal) = nan;
        pSelect(j,:,:) = p;
    end
    
    p_avg(numNodes==numNodesListFullyRec) = nanmean(pSelect(:));
    p_sd(numNodes==numNodesListFullyRec) = nanstd(pSelect(:)); 
    
    pRec(numNodes==numNodesListFullyRec) = mean(pRec_tmp); 
    pRecRandom(numNodes==numNodesListFullyRec) = ...
            p_avg(numNodes==numNodesListFullyRec)^numel(idNonDiagonal);     

    devFullyRec(numNodes==numNodesListFullyRec) = ...
                    pRec(numNodes==numNodesListFullyRec)./ ...
                    pRecRandom(numNodes==numNodesListFullyRec);
    
    fprintf('%d,%.2e,%.2e,%.2e\n',numNodes, ...
            pRec(numNodes==numNodesListFullyRec), ...
            pRecRandom(numNodes==numNodesListFullyRec), ...
            devFullyRec(numNodes==numNodesListFullyRec));
end

%% SAVE DATA
save([outputPath 'FullyRecurrentMotif_NumTrials_' num2str(numTrials) ...
    '_' postfix '.mat'], 'pRec','pRecRandom','numNodesListFullyRec', ...
    'p_avg','p_sd');