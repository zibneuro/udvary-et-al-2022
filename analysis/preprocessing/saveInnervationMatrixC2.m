% Read in InnervationMatrix csv files provided by Philipp Harth
% Saves Innervation as proper matrix incl. pre and postIDs
% for each cell type combination within C2 column
% in CellMatrix_C2
%   I.I: InnervationMatrix [Npre x Npost]
%   I.PreCellID [Npre x 1]
%   I.PostCellID [1 x Npost]
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

modelID = 'RBC_20';
matlabPath = 'D:\udvary-et-al-2022\analysis\';
load([matlabPath 'data\CellTable_' modelID '.mat'])
outputPath = [matlabPath 'data\CellMatrix_C2\'];
inputPath = [matlabPath 'preprocessing\data\DSC\'];

idx = strcmp(CellTable.nearest_column,'C2');
CellTable(~idx,:) = [];
[CellTypeList,~,ic] = unique(CellTable.cell_type); 
CellTypeCounts = accumarray(ic,1);
neuronIDList = cell(size(CellTypeList));

for ct = 1:length(CellTypeList)
    idxCT = strcmp(CellTable.cell_type,CellTypeList{ct});
    neuronIDList{ct} = CellTable.neuronID(idxCT);
end

clear idx idxCT ic;

%%
for preCT = 1:length(CellTypeList)-1
    
    precelltype = CellTypeList{preCT};
    preID_list = neuronIDList{preCT}';
    fprintf('Processing %s (n=%d) ...',precelltype,numel(preID_list));
    
    tbl = table; 
    for i = preID_list
        tbltmp = readtable([inputPath '\' num2str(i) '_DSC.csv'], ...
                    'Format','%d%f'); 
        tbltmp.pre_id = i.*int32(ones(size(tbltmp.post_id))); 
        tbl = [tbl; tbltmp]; 
    end
    clear tbltmp; 
    fprintf(' loaded all data ... ');

    %%
    for postCT = 1:length(CellTypeList)

        if strcmp(CellTypeList(postCT),'VPM')
            continue;
        end
        
        % Create Matrix
        I.PostCellID = int32(neuronIDList{postCT})';
        I.PreCellID = int32(neuronIDList{preCT});
        I.I = zeros(length(I.PreCellID),length(I.PostCellID)); 

        tmpID = ismember(tbl.post_id,neuronIDList{postCT}); 
        tbltmp = tbl(tmpID,:);
        
        [~,idxPre] = ismember(tbltmp.pre_id,I.PreCellID);
        [~,idxPost] = ismember(tbltmp.post_id,I.PostCellID);
        idx = sub2ind(size(I.I),idxPre,idxPost);
        I.I(idx) = tbltmp.DSC; 
        save([outputPath precelltype '_' CellTypeList{postCT} '.mat'],'I');
    end
    
    fprintf('saved matrices.\n');
end