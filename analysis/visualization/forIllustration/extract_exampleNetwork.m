% Extract example network of 50 nodes that has the 3 example cells in it
% only for C2 column
% only for illustration purposes!
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc
matlabPath = 'D:\udvary-et-al-2022\analysis\';
dataPath = [matlabPath 'data\example_selected_cells\'];
load([matlabPath 'data\CellTable_RBC_20.mat']);
idx = strcmp(CellTable.nearest_column,'C2');  
CellTable(~idx,:) = []; 

%%
L5PT_ID = 301854;
L2PY_ID = 748854;
L6ACC_ID = 199678;
numNodes = 50-3;
ids = []; 

for idtmp = [L5PT_ID L2PY_ID L6ACC_ID]
    preTbl = readtable([dataPath num2str(idtmp) '\' num2str(idtmp) ...
                        '_DSC.csv']);
    postTbl = readtable([dataPath num2str(idtmp) '\' num2str(idtmp) ...
                        '_innervating-pre-neurons-dsc.csv']);    
    ids = [ids; preTbl.post_id; postTbl.pre_id];
end

ids = unique(ids); 
idxtmp = ismember(ids,CellTable.neuronID); 
ids = ids(idxtmp); % all cells in C2 are involved somehow!

rng(11042021);
idxRnd = randperm(numel(ids),numNodes); 
idxtmp = ismember(CellTable.neuronID,[L5PT_ID L2PY_ID L6ACC_ID ids(idxRnd)']); 

exampleIDs = CellTable.neuronID(idxtmp);
exampleCTs = CellTable.cell_type(idxtmp); 

%% Create respective connectivity matrix
Imatrix = zeros(numel(exampleIDs)); 
CTlist = unique(exampleCTs);

for preCT = 1:length(CTlist)
    for postCT = 1:length(CTlist)
        load([matlabPath 'data\CellMatrix_C2\' CTlist{preCT} '_' ...
                        CTlist{postCT} '.mat']);
                    
        for preID = 1:length(exampleIDs)
             
            idxtmp_pre = exampleIDs(preID)==I.PreCellID;
            if sum(idxtmp_pre)==0
               continue; 
            end

            for postID = 1:length(exampleIDs)
                if preID==postID
                    continue;
                end
                
                idxtmp_post = exampleIDs(postID)==I.PostCellID;
                if sum(idxtmp_post)==0
                    continue;
                end
                
                Imatrix(preID,postID) = I.I(idxtmp_pre,idxtmp_post); 
            end
        end
        fprintf('%s-%s done\n',CTlist{preCT},CTlist{postCT});
    end  
end

pMatrix = 1-exp(-Imatrix); 

save([matlabPath 'data\NetworkGraphExample.mat'], ...
            'pMatrix','exampleIDs','exampleCTs');
