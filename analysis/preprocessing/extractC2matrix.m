% Saves .mat file of cell-to-cell DSC between all neurons in C2
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
load([matlabPath 'data\CellTable_RBC_20.mat'],'CellTable');
inputPath = [matlabPath 'preprocessing\data\'];
outputPath = [matlabPath 'preprocessing\data\'];

%% 
idx = strcmp(CellTable.nearest_column,'C2') & ...
            ~strcmp(CellTable.cell_type,'VPM'); 
neuronIDs = CellTable.neuronID(idx)';
clear CellTable;

% Create matrix
I = zeros(numel(neuronIDs));

%% 
for preID = 1:length(neuronIDs)
    tbl = readtable([inputPath 'DSC\' num2str(neuronIDs(preID)) '_DSC.csv']);
    for postID = 1:length(neuronIDs)
        idxtmp = (neuronIDs(postID)==tbl.post_id);
        if sum(idxtmp)==1
            I(preID,postID) = tbl.DSC(idxtmp); 
        end
    end
    
    fprintf('%d / %d\n',preID,length(neuronIDs));
end

%%
idx = logical(eye(size(I)));
I(idx) = nan; 
save([outputPath 'C2.mat'],'I','neuronIDs','-v7.3'); 