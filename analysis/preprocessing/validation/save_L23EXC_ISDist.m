%% Avermann et al.,2012
% Fig 7A
% L23->L23 Connection probability drops with increasing Intersomatic
% distance by -0.05 % per um
%
% Requires slices as .mat files (output from process_sliced_L23EXC.m)
% Generates: L23EXC_L23EXC_slice.mat
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
outputPath = [matlabPath 'data\'];
inputPath = [matlabPath 'preprocessing\data\L23EXC_L23EXC_slice\']; 
load([matlabPath 'data\CellTable_RBC_20.mat'],'CellTable');
idx = strcmp(CellTable.nearest_column,'C2');
neuronID = CellTable.neuronID(idx);
somaPos = [CellTable.soma_x(idx) CellTable.soma_y(idx) ...
            CellTable.soma_z(idx)];
clear CellTable;

%%
% Compute distance between soma positions
dim = [2]; % along y axis
ISLimit = [0:40:160];
ISLimit_center = ISLimit(1:end-1)+(ISLimit(2)-ISLimit(1))/2;

CP_median = nan(1,length(ISLimit)-1);
CP_25 = nan(size(CP_median));
CP_75 = nan(size(CP_median));

for i = 1:length(ISLimit)-1
                    
    ISLimit_y = [ISLimit(i) ISLimit(i+1)];
    p = [];

    for sliceID = 0:19 
        load([inputPath 'slice-' num2str(sliceID) '.mat'],'I');
        distMatrix = nan(size(I.I)); 
        n = size(I.I,1); 
        
        for c1 = 1:n
            idx1 = ismember(neuronID,I.PreCellID(c1));
            for c2 = 1:c1
                idx2 = ismember(neuronID,I.PreCellID(c2));
                distMatrix(c1,c2) = sqrt(sum((somaPos(idx1,dim) ...
                                            - somaPos(idx2,dim)).^2));
                distMatrix(c2,c1) = distMatrix(c1,c2);
            end
        end
        
        idx = distMatrix>ISLimit_y(1) & distMatrix<=ISLimit_y(2);
        p = [p; 1-exp(-I.I(idx))]; 
        fprintf('%d (%.2f)\n',sliceID,mean(ISLimit_y));
    end
                                           
    CP_median(i) = median(p);
    CP_25(i) = prctile(p,25);
    CP_75(i) = prctile(p,75);
    
    fprintf('%.2f done\n',mean(ISLimit_y));
end 

%%
save([outputPath '\L23EXC_L23EXC_slice.mat'], ...
            'CP_median','CP_25','CP_75','ISLimit','ISLimit_center');
