%% Display number of Cells within Layer 5
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc
matlabPath = 'D:\udvary-et-al-2022\analysis\';

%% Ratio of cells of each cell type in C2 in L5
load([matlabPath 'data\CellTable_RBC_20.mat']);
idxtmp = strcmp(CellTable.nearest_column,'C2');
CellTable(~idxtmp,:) = [];

CellTypeList = unique(CellTable.cell_type);

fprintf('Cells of that Cell Type in Layer 5\n');
for i = 1:length(CellTypeList)
   idx_ct = strcmp(CellTable.cell_type,CellTypeList{i});
   idx_L5 = (CellTable.layer(idx_ct)==5);
   
   if sum(idx_L5)>0
        fprintf('%s: %d/%d (%.2f%%)\n',CellTypeList{i},sum(idx_L5), ...
            sum(idx_ct),sum(idx_L5)/sum(idx_ct)*100);
   end
end

% Excitatory cells in Layer 5
fprintf('\nCell type composition of Layer 5\n');
nL5 = sum(CellTable.layer==5 & ~strcmp(CellTable.cell_type,'INH'));
for i = 2:length(CellTypeList)
   idx_ct = strcmp(CellTable.cell_type,CellTypeList{i});
   idx_L5 = (CellTable.layer(idx_ct)==5);
   
   if sum(idx_L5)>0
        fprintf('%s: %d/%d (%.2f%%)\n',CellTypeList{i},sum(idx_L5), ...
            sum(idx_ct),sum(idx_L5)/nL5*100);
   end
end