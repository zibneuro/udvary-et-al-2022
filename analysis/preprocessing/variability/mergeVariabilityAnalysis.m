%% Merge Variability Analysis
% Output: MergedAxonDend.mat
% Requires: Axon.mat and Dendrite.mat (from
% variabilityAnalysis_DiffMorphSample.m)
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
inputPath = [matlabPath 'preprocessing\data\variability\']; 
outputPath = [matlabPath 'data\variability\']; 

numVoxels = 512; 
axon = load([inputPath 'Axon.mat'],'results'); 
dend = load([inputPath 'Dendrite.mat'],'results'); 

%%
CV_Len = nan(1,length(axon.results)-1);
CV_CC = nan(1,length(axon.results)-1);

for i = 1:length(axon.results)-1

    % NumVoxels x Combinations
    axon_lentmp = axon.results{i}.len;
    axon_cctmp = axon.results{i}.contributingCells;      
    dend_lentmp = dend.results{i}.len;
    dend_cctmp = dend.results{i}.contributingCells;  
    
    Ncombs = size(axon_lentmp,2);
    
    if size(axon_lentmp,2)~=size(dend_lentmp,2)
       error('Not same sample size (length)'); 
    end
    if size(axon_cctmp,2)~=size(dend_cctmp,2)
       error('Not same sample size (cc)'); 
    end
    if size(axon_cctmp,2)~=size(axon_lentmp,2)
       error('Not same sample size (len and cc)'); 
    end
        
    % All combinations between axon and dendrite combinations
    % -> 500 x 500
    % -> sort by cube
    % -> NumVoxels x Combinations (500x500)
    lentmp = nan(numVoxels,Ncombs*Ncombs); 
    for v = 1:numVoxels
        x = repmat(axon_lentmp(v,:),Ncombs,1) + ...
            repmat(dend_lentmp(v,:)',1,Ncombs);
        lentmp(v,:) = x(:); 
    end
    
    cctmp = nan(numVoxels,Ncombs*Ncombs); 
    for v = 1:numVoxels
        x = repmat(axon_cctmp(v,:),Ncombs,1) + ...
            repmat(dend_cctmp(v,:)',1,Ncombs);
        cctmp(v,:) = x(:); 
    end
    
    % Over combinations -> NumVoxels x 1
    % Only save median over combinations
    tmp = std(lentmp,[],2)./mean(lentmp,2);
    CV_Len(i) = median(tmp);
    
    tmp = std(cctmp,[],2)./mean(cctmp,2);
    CV_CC(i) = median(tmp);
    
    fprintf('%d done\n',i);
end

%%
save([outputPath 'MergedAxonDend.mat'],'CV_CC','CV_Len'); 