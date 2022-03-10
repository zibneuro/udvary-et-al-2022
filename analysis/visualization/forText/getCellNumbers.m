%% Display Cell Numbers (composition of model)
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
load([matlabPath 'data\CellTable_RBC_20.mat'],'CellTable');

idxC2 = strcmp(CellTable.nearest_column,'C2'); 
idxVPM = strcmp(CellTable.cell_type,'VPM');
idxINH = strcmp(CellTable.cell_type,'INH');
idxEXC = ~idxVPM & ~idxINH;

fprintf('VPM = %d\nINH = %d\nEXC = %d\n',sum(idxVPM),sum(idxINH),sum(idxEXC));
fprintf('VPM[C2] = %d\nINH[C2] = %d\nEXC[C2] = %d\n',sum(idxVPM & idxC2), ...
            sum(idxINH & idxC2),sum(idxEXC & idxC2));
        
% Layer 5 stats
idxL5 = (CellTable.layer==5);
[celltypes_in_l5,~,ic] = unique(CellTable.cell_type(idxL5 & idxC2));
a_counts = accumarray(ic,1);

for i = 1:length(celltypes_in_l5)
    fprintf('%s[C2,L5] = %d (%.2f%%)\n',celltypes_in_l5(i),a_counts(i), ...
        a_counts(i)/sum(a_counts)*100);
end

