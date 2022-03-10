%% Perin et al., 2011
% Fig 1E
% L5tt -> L5tt Connection Probability drops with increasing Intersomatic 
% distance
% ISLimit_center = ([66.557 77.846 88.665 99.601 110.537 121.708 132.88 ...
%                             143.933 154.752]-61.009)./(155.692-61.009).*300;
% rangePerin = 25; 
%
% Requires slices as .mat files (output from process_sliced_L5PT.m)
% Generates: L5PT_L5PT_slice.mat
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
outputPath = [matlabPath 'data\'];
inputPath = [matlabPath 'preprocessing\data\L5PT_L5PT_slice\']; 
load([matlabPath 'data\CellTable_RBC_20.mat'],'CellTable');
idx = strcmp(CellTable.nearest_column,'C2') & ...
            strcmp(CellTable.cell_type,'L5PT');

neuronID = CellTable.neuronID(idx);
somaPos = [CellTable.soma_x(idx) CellTable.soma_y(idx) ...
            CellTable.soma_z(idx)];
clear CellTable;

%%
% Compute distance between soma positions
dim = [2]; % along y axis
ISLimit = [0:35:330];
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
        
    end
                                           
    CP_median(i) = median(p);
    CP_25(i) = prctile(p,25);
    CP_75(i) = prctile(p,75);
    
    fprintf('%.2f done\n',mean(ISLimit_y));
end 

%%
save([outputPath 'L5PT_L5PT_slice.mat'], ...
            'CP_median','CP_25','CP_75','ISLimit','ISLimit_center');
