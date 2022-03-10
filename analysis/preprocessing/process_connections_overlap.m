%% Calculation for
% Occurrences of neuron pairs vs. connected overlap volumes
% Outputs: dsc_overlapping-cubes_batch-*.mat files
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
outputPath = [matlabPath 'data\dataPetersRuleBatch\']; 
inputPath = [matlabPath 'preprocessing\data\connections_overlap\']; 

%%
edge_sz = 0.01;
perc_edges = [0:edge_sz:1];
batchID_list = [0:9];
numBins = length(perc_edges)+1;
% 1 -> 0%
% N -> 100%
% >1 <N other bins

%%
%parpool(6);
%parfor batchID = batchID_list
for batchID = batchID_list(1)
    tbl = readtable([inputPath 'dsc_overlapping-cubes_batch-' ...
                            num2str(batchID)]);
    % tbl.Var1: dsc
    % tbl.Var2: number of overlapping cubes
    Npairs = size(tbl,1); 
        
    perc_p_mean = zeros(1,numBins);
    perc_p_sd = zeros(1,numBins);
    perc_p_min = 2.*ones(1,numBins);
    perc_p_max = -ones(1,numBins);
    perc_p_N = zeros(1,numBins);
    
    % No synapses across overlap volumes // Bin 1
    idx = 1; 
    xtmp = poisspdf(0,tbl.Var1);
    perc_p_mean(idx) = mean(xtmp);
    perc_p_sd(idx) = std(xtmp);
    perc_p_min(idx) = min(xtmp);
    perc_p_max(idx) = max(xtmp);
    perc_p_N(idx) = numel(xtmp);
    
    % As many synapses as overlap volumes // Bin N
    idx = numBins; 
    xtmp = poisspdf(tbl.Var2,tbl.Var1);
    perc_p_mean(idx) = mean(xtmp);
    perc_p_sd(idx) = std(xtmp);
    perc_p_min(idx) = min(xtmp);
    perc_p_max(idx) = max(xtmp);
    perc_p_N(idx) = numel(xtmp);
    
    % All in between 0% and 100%
    maxOverlap = max(tbl.Var2)-1;
    for i = 1:maxOverlap
        idxDel = (tbl.Var2==i);
        tbl(idxDel,:) = []; 
        xtmp = poisspdf(i,tbl.Var1);
        idx = ceil(i./tbl.Var2./edge_sz)+1;
        % ceil([0.005 0.01 0.02 0.999]./0.01)+1
        %       -> [2 2 3 101]
        if i<maxOverlap/4
            fprintf('%d / %d (n=%d) (batchID = %d)\n', ...
                        i,maxOverlap,numel(xtmp),batchID);
        end
        idxList = unique(idx)';
        for j = idxList
            idxtmp = (j==idx);
            xtmp2 = xtmp(idxtmp);
            Ntmp2 = sum(idxtmp);
            sdtmp2 = std(xtmp2);
            mtmp2 = mean(xtmp2);
            
            if perc_p_N(j)>0
                [M,SD] = totalMSD([perc_p_mean(j) mtmp2], ...
                                [perc_p_sd(j) sdtmp2], ...
                                [perc_p_N(j) Ntmp2]);
                perc_p_sd(j) = SD;
                perc_p_mean(j) = M;
            else
                perc_p_sd(j) = sdtmp2;
                perc_p_mean(j) = mtmp2;
            end
            perc_p_min(j) = min([perc_p_min(j); xtmp2]);
            perc_p_max(j) = max([perc_p_max(j); xtmp2]);
            perc_p_N(j) = perc_p_N(j) + Ntmp2;
        end
    end    

    % Check data
    if sum(perc_p_min>1)>0 || sum(perc_p_max<0)>0 || sum(perc_p_N==0)>0
       error('something went wrong in batch %d',batchID); 
    end   
    
    fprintf('saving batchID %d\n',batchID);

    result = struct; 
    result.perc_p_mean = perc_p_mean;
    result.perc_p_sd = perc_p_sd;
    result.perc_p_min = perc_p_min;
    result.perc_p_max = perc_p_max;
    result.perc_p_N = perc_p_N;
    result.edge_sz = edge_sz;
    result.perc_edges = perc_edges; 
    result.Npairs = Npairs; 
    
    parsave([outputPath 'dsc_overlapping-cubes_batch-' ...
                num2str(batchID) '.mat'],result);
end
     
function parsave(fname,result)
    save(fname,'result');
end