%% Calculates doublet motif distribution across L5PT slices
% Requires slices as .mat files (output from process_sliced_L5PT.m)
% Outputs: L5PT_Doublets.mat
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
outputPath = [matlabPath 'data\'];
inputPath = [matlabPath 'preprocessing\data\L5PT_L5PT_slice\'];

%%
sliceIDs = [0:1:19];
doubletMotifLabel = {'uni','bi','unconn'}; 
doubletMotifP_RND = nan(length(sliceIDs),3);
doubletMotifP = nan(length(sliceIDs),3);
numSamples = nan(1,length(sliceIDs)); 

for k = 1:length(sliceIDs)
   
    load([inputPath 'slice-' num2str(sliceIDs(k)) '.mat']);
    numCells = numel(I.PreCellID);
    N = numCells*numCells-numCells;

    p = 1-exp(-I.I); 
    ptmp = p;
    ptmp(logical(eye(size(ptmp)))) = nan; % ignore diagonal
    p_avg = nanmean(ptmp(:));

    % Two Neuron Connectivity Patterns
    % p is probability of being connected
    %   unidirectionally connected: p*(1-p) + (1-p)*p
    %   bidirectionally connected: p^2
    %   unconnected: (1-p)^2
    % unidirectional, bidirectional, unconnected
    doubletMotifP_RND(k,:) = [2*p_avg*(1-p_avg) p_avg^2 (1-p_avg)^2];
    
    % Sample
    twoMotifs = nan(N,3);
    c = 1; 
    for i = 1:numCells
        for j = 1:numCells
            if i==j
               continue; 
            end
            p1 = p(i,j);
            p2 = p(j,i);
            twoMotifs(c,:) = [p1*(1-p2)+p2*(1-p1) p1*p2 (1-p1)*(1-p2)]; 
            c = c+1; 
        end
    end
    
    doubletMotifP(k,:) = mean(twoMotifs);
    numSamples(k) = N; 
end

doubleMotifDev = doubletMotifP./doubletMotifP_RND; 
save([outputPath 'L5PT_Doublets.mat'], ...
        'doubleMotifDev','doubletMotifP','doubletMotifP_RND', ...
        'sliceIDs','doubletMotifLabel','numSamples');