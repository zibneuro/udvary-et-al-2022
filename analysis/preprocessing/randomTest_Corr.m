%% Perform random permuation test for correlation between predicted and 
% measured connection probabilities
% Outputs: correlationValuesRandom.mat
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
outputPath = [matlabPath 'data\']; 

%%
tbl = readtable([matlabPath 'data\cellular_connectivity\stats\summary.csv']); 
CP_measured = tbl.empirical';
CP_AVG_predicted = tbl.avg_slices_merged';
[r,p] = corr(CP_AVG_predicted',CP_measured');

%% Rand Perm Test
rng(1);
numTrails = 1e5;
rRandom = nan(1,numTrails);
rej = 0; 

for t = 1:numTrails
    idxRnd = randperm(length(CP_AVG_predicted));
    rRandom(t) = corr(CP_AVG_predicted(idxRnd)',CP_measured');
    
    if rRandom(t)>r
        rej = rej + 1; 
    end
end

save([outputPath 'correlationValuesRandom.mat'],'r','rRandom','numTrails','rej');