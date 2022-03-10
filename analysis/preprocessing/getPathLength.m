%% Extracts path length values across subvolumes
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc
warning on;

matlabPath = 'D:\udvary-et-al-2022\analysis\';
load([matlabPath 'data\CellTable_RBC_20.mat'],'CellTable');
dataPath1 = [matlabPath 'data\']; 
dataPath2 = [matlabPath 'preprocessing\data\subcellular_features\'];
outputPath = [matlabPath 'data\dataSubvolume\']; 

% Delete VPM
idxVPM = strcmp(CellTable.cell_type,'VPM');
CellTable(idxVPM,:) = []; 

% CHANGE subvolume dimensions here!
%volStr = '100-100-50';
volStr = '50-50-50';

%% Get grid
tbl_grid = readGrid([dataPath1 'grid_' volStr '_ref-volume.csv']);
pos = [tbl_grid.ix tbl_grid.iy tbl_grid.iz];

%% Do for axon and dendrite
for k = 1:2
    
    if k==1
        synSide = 'post';
    else
        synSide = 'pre';
    end
    
    % Get list of innervating CellIDs
    filename = ['innervating_ref-volume_' synSide '.txt']; 
    fileID = fopen([dataPath2 filename],'r');
    ids = cell2mat(textscan(fileID,'%u32'));
    idxtmp = ismember(ids,CellTable.neuronID);
    ids = ids(idxtmp);
    fprintf('Removed %d cells outside of vS1 (%d/%d)\n', ...
            numel(idxtmp)-sum(idxtmp),sum(idxtmp),numel(idxtmp));
    
    % iterate over all innervating cells
    len = [];
    tmp = round(length(ids)/500);
    for i = 1:length(ids)
        tbl_tmp = readtable([dataPath2 ...
                'subcellular_features_' synSide 'synaptic_' volStr '\' ...
                num2str(ids(i)) '.csv']);
        postmp = [tbl_tmp.ix tbl_tmp.iy tbl_tmp.iz];
        idxtmp = ismember(postmp,pos,'rows');

        if sum(idxtmp)>0
            len = [len; tbl_tmp.distance_soma(idxtmp)]; 
        else
           warning('No overlap in volume found! ID = %d',ids(i)); 
        end

        if rem(i,tmp)==0
            fprintf('%.1f%% (%s)\n',i/length(ids)*100,synSide);
        end
    end

    % Save pathlength
    save([outputPath synSide '_pathLength_' strrep(volStr,'-','_') ...
                '.mat'],'len');
end