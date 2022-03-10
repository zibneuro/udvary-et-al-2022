% Generates *_Axon.mat and *_Dendrite.mat for each soma distribution
% 4 soma distributions (RBC_21,22,23,24)
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
dataPath = [matlabPath 'data\']; 
outputPath = [matlabPath 'preprocessing\data\variability\']; 
tblGrid = readtable([dataPath '\grid_50-50-50_ref-volume.csv']); 
numVoxels = size(tblGrid,1); 
idList = 21:24;               

%% Iterate over post and pre
for id = idList
    for k = 1:2
        if k==1
            str = 'post';
            strtype = 'Dendrite';
            listMorph = 1:30;
            listReal = 0:499;
        else
            str = 'pre';
            strtype = 'Axon';
            listMorph = 1:14;
            listReal = 0:499;
        end

        filepath = [outputPath 'RBC_' num2str(id) '\' str '\'];
        results = cell(1,max(listMorph)); 

        for i = listMorph

            % Features [numVoxels x numReal]
            features.len = zeros(numVoxels,length(listReal));
            features.contributingCells = zeros(numVoxels,length(listReal));
            features.boutons = zeros(numVoxels,length(listReal));
            features.branches = zeros(numVoxels,length(listReal));
            features.path_soma = zeros(numVoxels,length(listReal));
            features.cell_types = zeros(numVoxels,length(listReal));

            for j = 1:length(listReal)
                filename = [filepath str '_morphologies-' num2str(i) ...
                    '_realization-' num2str(listReal(j))]; 
                tbl = readtable(filename);

                for v = 1:numVoxels
                    idxtmp = (tblGrid.ix(v)==tbl.ix & tblGrid.iy(v)==tbl.iy & ...
                                tblGrid.iz(v)==tbl.iz); 

                    if sum(idxtmp)~=1
                        error('volume not found?!');
                    end

                    features.len(v,j) = tbl.length(idxtmp); 
                    features.contributingCells(v,j) = ...
                                tbl.contributing_cells(idxtmp); 
                    features.boutons(v,j) = tbl.boutons(idxtmp);
                    features.branches(v,j) = tbl.branches(idxtmp);
                    features.path_soma(v,j) = tbl.path_soma_sum(idxtmp);
                    features.cell_types(v,j) = tbl.cell_types(idxtmp);
                end
            end

            results{i} = features; 
            fprintf('%d_%s: %d/%d done\n',id,strtype,i,max(listMorph));
        end

        save([outputPath 'RBC_' num2str(id) '_' strtype '.mat'],'results');
    end
end