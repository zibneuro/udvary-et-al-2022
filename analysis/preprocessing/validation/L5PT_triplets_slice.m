%% Calculates triplet motif distribution across L5PT slices
% Requires slices as .mat files (output from process_sliced_L5PT.m)
% Outputs: L5PT_Triplets.mat
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
outputPath = [matlabPath 'data\'];
inputPath = [matlabPath 'preprocessing\data\L5PT_L5PT_slice\'];

MotifSpectrum = getMotifSpectrum();
numMotifs = numel(MotifSpectrum.Class); 

%%
sliceIDs = [0:1:19];
numTrials = 1e4; % 10,000 per slice
rng(50000); 
comb = [1 2; 1 3; 2 1; 2 3; 3 1; 3 2]'; 
pMotif64 = nan(length(sliceIDs),64);
p_avg = zeros(length(sliceIDs),1);

for k = 1:length(sliceIDs)
   
    load([inputPath 'slice-' num2str(sliceIDs(k)) '.mat']);
    p = 1-exp(-I.I); 
    pMeanMotif64_tmp = zeros(1,numMotifs); 
    N = length(I.PreCellID); 

    prevProgress = 0; 
    fprintf('>> 0%%\n');
    
    for t = 1:numTrials
        idx_tmp = randperm(N,3);
        ptmp = zeros(3,3); 

        for i = 1:size(comb,2)
           idxMatrix = comb(:,i);
           idx1 = idx_tmp(idxMatrix);
           ptmp(idxMatrix(1),idxMatrix(2)) = p(idx1(1),idx1(2)); 
        end

        ptmp(logical(eye(3,3))) = nan; 
        p_avg(k) = p_avg(k) + nanmean(ptmp(:)); 

        for motifIDX = 1:numMotifs
            m = squeeze(MotifSpectrum.Connectivity(motifIDX,:,:));
            m(logical(eye(3,3))) = nan; 
            pMeanMotif64_tmp(motifIDX) = pMeanMotif64_tmp(motifIDX) + ...
                                     prod(ptmp(m==1)) * prod(1-ptmp(m==0)); 
        end

        % Display progress
        progress = floor(t/numTrials*100);
        if (progress~=prevProgress || t==numTrials)
            if (floor(t/numTrials*100))>9
                str = '\b\b\b';
            else
                str = '\b\b';
            end
            fprintf([str '%d%%'],progress);  
        end
        prevProgress = progress;

    end
    fprintf(' Slice %d done\n',sliceIDs(k));
    pMotif64(k,:) = pMeanMotif64_tmp./numTrials; 
    p_avg(k) = p_avg(k)./numTrials;
    
end

save([outputPath 'L5PT_Triplets.mat'],'p_avg','pMotif64','numTrials');