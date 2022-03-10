%% Compute how many axons/dendrites originate from cells whose somata is 
% inside cube, 
% outside cube (but within vS1), 
% outside cube (but from VPM),
% outside cube (outside of vS1 and not in VPM),
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
dataPath1 = [matlabPath 'data\'];
load([matlabPath 'data\CellTable_RBC_20.mat'],'CellTable');
dataPath2 = [matlabPath 'preprocessing\data\subcellular_features\'];
outputPath = [matlabPath 'data\dataSubvolume\']; 

% CHANGE subvolume dimensions here!
%volStr = '100-100-50';
%cubesz = [100 100 50];
volStr = '50-50-50'; 
cubesz = [50 50 50];

tbl_grid = readGrid([dataPath1 'grid_' volStr '_ref-volume.csv']);
cubeID = [tbl_grid.ix tbl_grid.iy tbl_grid.iz];
minBB = [tbl_grid.cube_origin_x tbl_grid.cube_origin_y ...
            tbl_grid.cube_origin_z];
maxBB = minBB + cubesz;
numCubes = size(cubeID,1);

%%
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
    
    % iterate over all innervating cells
    % Inside,Outside(vS1),Outside(VPM),Outside
    cellInsideOut = zeros(numCubes,4); 
    tmp = round(length(ids)/500);
    for i = 1:length(ids)
        
        tbl_tmp = readtable([dataPath2 ...
                'subcellular_features_' synSide 'synaptic_' volStr '\' ...
                num2str(ids(i)) '.csv']);            
        cubeIDtmp = [tbl_tmp.ix tbl_tmp.iy tbl_tmp.iz];
        cube_bool = ismember(cubeID,cubeIDtmp,'rows');
        
        if sum(cube_bool)==0
           warning('No overlap in volume found! ID = %d',ids(i)); 
           continue;
        end
        
        idxNeuron = (ids(i)==CellTable.neuronID);
        if sum(idxNeuron)==1

            % in VPM -> outside of cube
            if strcmp(CellTable.cell_type(idxNeuron),'VPM')
                idxtmp = 3;
                cellInsideOut(cube_bool,idxtmp) = ...
                        cellInsideOut(cube_bool,idxtmp)+1; 
            else % other cell types (iterate over all cubes)
                somaPos = [CellTable.soma_x(idxNeuron) ...
                        CellTable.soma_y(idxNeuron) ...
                        CellTable.soma_z(idxNeuron)];
                
                 cubeIDX = find(cube_bool)';
                 for ii = cubeIDX
                     % Check whether soma is within cube or outside
                     if sum(minBB(ii,:)<=somaPos & maxBB(ii,:)>=somaPos)==3
                        cellInsideOut(ii,1) = cellInsideOut(ii,1)+1; 
                     else
                        cellInsideOut(ii,2) = cellInsideOut(ii,2)+1; 
                     end
                 end 
            end
            
        else % outside of vS1+VPM
            idxtmp = 4;
            cellInsideOut(cube_bool,idxtmp) = ...
                    cellInsideOut(cube_bool,idxtmp)+1; 
        end

        if rem(i,tmp)==0
            fprintf('%.1f%% (%s)\n',i/length(ids)*100,synSide);
        end
    end
    
    % Save pathlength
    save([outputPath synSide ...
            '_CellsInsideOut_' strrep(volStr,'-','_') '.mat'], ...
            'cellInsideOut');
end