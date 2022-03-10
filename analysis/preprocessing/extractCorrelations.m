% Go through all cell type combiations and calculates their inDegree correlations
% Outputs: CorrData.mat
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
dataPath = [matlabPath 'data\']; 
load([matlabPath 'data\CellTable_RBC_20.mat'],'CellTable');
idx = strcmp(CellTable.nearest_column,'C2');
CellTypeList = unique(CellTable.cell_type); 
clear CellTable;

%%   
PreType1 = cell(0);
PreType2 = cell(0);
PostType = cell(0);
R = []; 
P = []; 

for idxCTPre1 = 1:length(CellTypeList)

    for idxCTPre2 = 1:length(CellTypeList)
        
        if idxCTPre1<=idxCTPre2
            continue;
        end
                
        for idxCTPost = 1:length(CellTypeList)
                        
            if strcmp(CellTypeList{idxCTPost},'VPM') % Skip VPM
                continue;
            end
            
            load([dataPath 'CellMatrix_C2\' ...
                CellTypeList{idxCTPre1} '_' ...
                CellTypeList{idxCTPost} '.mat'],'I');
            I1 = I;
            load([dataPath 'CellMatrix_C2\' ...
                CellTypeList{idxCTPre2} '_' ...
                CellTypeList{idxCTPost} '.mat'],'I');
            I2 = I;
        
            if sum(I1.PostCellID==I2.PostCellID)~=numel(I2.PostCellID)
               error('PostCellIDs not sorted!'); 
            end
            
            syn1 = sum(I1.I,1);
            syn2 = sum(I2.I,1); 
            [r,p] = corr(syn1',syn2');            
            R(end+1) = r; 
            P(end+1) = p;
            PreType1{end+1} = CellTypeList{idxCTPre1};
            PreType2{end+1} = CellTypeList{idxCTPre2};
            PostType{end+1} = CellTypeList{idxCTPost}; 
        end   
        fprintf('%s-%s done\n',CellTypeList{idxCTPre1}, ...
                            CellTypeList{idxCTPre2});
    end
end

save([dataPath 'CorrData.mat'],'R','P','PreType1','PreType2','PostType');