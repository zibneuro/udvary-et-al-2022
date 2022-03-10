% Run motif analysis for chain motifs of 3 to 10 nodes
% Requires: C2.mat (generated from extractC2matrix.m)
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
inputPath = [matlabPath 'preprocessing/data/']; 
outputPath = [matlabPath 'data/higherOrder/']; 
load([inputPath 'C2.mat'],'I'); 
postfix = 'C2';
pTotal = 1-exp(-I);
clear I; 
numTrials = 10000;
numChainMotifs = 1000; 
numNodesList = [3:1:10]; 
n = size(pTotal,1); 

%%
rng(154581); 
for numNodes = numNodesList
    
    edgeMatrix = generateRandomChainMotif(numNodes,numChainMotifs);
    pSelect = nan(numTrials,numNodes,numNodes); 
    idDiagonal = 1:numNodes+1:(numNodes^2); % id of diagonal values

    for j = 1:numTrials
        idx = randperm(n,numNodes);
        p = pTotal(idx,idx);
        p(idDiagonal) = nan; 
        pSelect(j,:,:) = p; 
    end

    p_avg = nanmean(pSelect(:));
    p_sd = nanstd(pSelect(:)); 
    
    [edgeStats] = extractEdgeMotifProbabilityGivenChainMotif( ...
                                edgeMatrix,pSelect,p_avg);     
    edgeStats.comment = 'chain'; 
    edgeStats.p_sd = p_sd; 

    % SAVE DATA
    save([outputPath 'ChainMotif_NumTrials_' num2str(numTrials) ...
        '_NumCells_' num2str(numNodes) '_' postfix '.mat'], ...
        'edgeStats');
    fprintf('SAVED for %d\n',numNodes); 
end