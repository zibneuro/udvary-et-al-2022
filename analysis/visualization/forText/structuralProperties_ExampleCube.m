%% compute structural composition of example cube (L2PY+L5PT)
% - number of branchlets in cube
% - number of contributing cells in cube
%
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
filepath1 = [matlabPath 'data\'];  
filepath2 = [matlabPath 'data\dataSubvolume\'];
tblGrid = readtable([filepath1 'grid_50-50-50_ref-volume.csv']); 
    
%%
% example cube
% [-100 300 350] [-50 350 400]
% -> [-75 325 375]
targetCube = [-75 325 375];
idxCube = tblGrid.cube_origin_x==-100 & tblGrid.cube_origin_y==300 & ...
                tblGrid.cube_origin_z==350;
            
filename = 'cube-stats_ref-volume\cube-stats_50-50-50_ref-volume_pre.csv';
tblAxon = readtable([filepath1 filename]);
idxAxon = tblAxon.ix==tblGrid.ix(idxCube) & tblAxon.iy==tblGrid.iy(idxCube) & ...
                tblAxon.iz==tblGrid.iz(idxCube);
tblAxon(~idxAxon,:) = [];

% Dendrite / Post
filename = 'cube-stats_ref-volume\cube-stats_50-50-50_ref-volume_post.csv';
tblDend = readtable([filepath1 filename]);
idxDend = tblDend.ix==tblGrid.ix(idxCube) & tblDend.iy==tblGrid.iy(idxCube) & ...
                tblDend.iz==tblGrid.iz(idxCube);
tblDend(~idxDend,:) = [];

%% Display branchlets in cube
fprintf('#DendBranchlets = %.2f\n#AxonBranchlets = %.2f\n', ...
    tblDend.branches,tblAxon.branches);
fprintf('#Boutons = %.0f\n',tblAxon.boutons);

%%
% Axon Contributing Cells (same sorting as tblGrid)
load([filepath2 'pre_CellTypeContribution_50_50_50.mat']);
numPreCells = sum(ctContribution(idxCube,:));

load([filepath2 'post_CellTypeContribution_50_50_50.mat']);
numPostCells = sum(ctContribution(idxCube,:));

fprintf('#PostCells = %.2f\n#PreCells = %.2f\n#TotelCells = %.2f\n', ...
    numPostCells,numPreCells,numPreCells+numPostCells);