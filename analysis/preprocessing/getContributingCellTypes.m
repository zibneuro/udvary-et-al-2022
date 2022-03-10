%% Extracts number of cell types contributing/innervating to each cube
% calculates for pre- and postsynaptic neurons (i.e., axons and dendrites)
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc
warning on;

matlabPath = 'D:\udvary-et-al-2022\analysis\';
dataPath1 = [matlabPath 'data\']; 
dataPath2 = [matlabPath 'preprocessing\data\subcellular_features\'];
outputPath = [matlabPath 'data\dataSubvolume\']; 

% CHANGE subvolume dimensions here!
%volStr = '100-100-50';
volStr = '50-50-50';

%% Get grid
tbl_grid = readGrid([dataPath1 'grid_' volStr '_ref-volume.csv']);
pos = [tbl_grid.ix tbl_grid.iy tbl_grid.iz];
numCubes = size(pos,1);

%% CellTypes
load([matlabPath 'data\CellTable_RBC_20.mat'],'CellTable');
% Exc Cell Types
% Inh Layer 1,2,3,4,5,6
idxINH = strcmp(CellTable.cell_type,'INH');
tmp = strcat('INH_L',num2str(CellTable.layer(idxINH)));
tmp = cellstr(tmp);
tmp = strrep(tmp,' ','');
CellTable.cell_type(idxINH) = tmp; 
CellTypeList = unique(CellTable.cell_type); 

%% Do for axon and dendrite
for k = 1:2
    
    if k==1
        synSide = 'post';        
    else
        synSide = 'pre';
    end
    skippedCells = 0; 
    
    % Get list of innervating CellIDs
    filename = ['innervating_ref-volume_' synSide '.txt']; 
    fileID = fopen([dataPath2 filename],'r');
    ids = cell2mat(textscan(fileID,'%u32'));

    % iterate over all innervating cells
    ctContribution = zeros(numCubes,length(CellTypeList));
    tmp = round(length(ids)/500);
    for i = 1:length(ids)
        
        idxtmp = (ids(i)==CellTable.neuronID);
        if sum(idxtmp)~=1
           %warning('CellID %d found %d in Table!',ids(i),sum(idxtmp)); 
           % Probably this cell is not within vS1
           skippedCells = skippedCells + 1; 
           continue;
        end
        idxCT = strcmp(CellTypeList,CellTable.cell_type(idxtmp));

        tbl_tmp = readtable([dataPath2 ...
                'subcellular_features_' synSide 'synaptic_' volStr '\' ...
                num2str(ids(i)) '.csv']);            
        postmp = [tbl_tmp.ix tbl_tmp.iy tbl_tmp.iz];
        idxCube = ismember(pos,postmp,'rows');
        
        if sum(idxCube)==0
           warning('No overlap in volume found! ID = %d',ids(i)); 
        end

        ctContribution(idxCube,idxCT) = ctContribution(idxCube,idxCT)+1; 

        if rem(i,tmp)==0
            fprintf('%.1f%% (%s)\n',i/length(ids)*100,synSide);
        end
    end
    
    fprintf('Skipped: %d / %d (%.2f%%)\n',skippedCells,length(ids), ...
        skippedCells/length(ids)*100);

    % Save pathlength
    save([outputPath synSide ...
            '_CellTypeContribution_' strrep(volStr,'-','_') '.mat'], ...
            'ctContribution','CellTypeList','pos');
end