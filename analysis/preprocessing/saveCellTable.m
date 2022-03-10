% Creates CellTable_RBC_20
% table contains meta information about neurons in model
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
modelID = 'RBC_20';
inputPath = [matlabPath 'preprocessing\data\' modelID '\'];
outputPath = [matlabPath 'data\'];
addpath(genpath([matlabPath 'functions\']));

%%
tblRegions = readRegionsCSV([inputPath 'regions.csv']);
tblCellTypes = readCellTypesCSV([inputPath 'cell_types.csv']);
CellTable = readNeuronsCSV([inputPath 'neurons.csv']);

CellTable.cell_type = tblCellTypes.name(CellTable.cell_type+1);
idxVPM = strcmp(CellTable.cell_type,'VPM');
idxDel = ~(CellTable.inside_vS1 | idxVPM); % Delete cells outside
CellTable(idxDel,:) = [];
CellTable.inside_vS1 = []; 
CellTable.synaptic_side = []; 
CellTable.laminar_location = []; 
CellTable.Properties.VariableNames(1) = {'neuronID'};
idxVPM = CellTable.nearest_column==-1; 
CellTable.nearest_column(idxVPM) = CellTable.region(idxVPM); 
CellTable.nearest_column = tblRegions.name(CellTable.nearest_column+1);
CellTable.region = tblRegions.name(CellTable.region+1);
CellTable.cortical_depth(idxVPM) = nan; 
CellTable.soma_x(idxVPM) = nan; 
CellTable.soma_y(idxVPM) = nan; 
CellTable.soma_z(idxVPM) = nan; 

% Delete barreloid
CellTable.nearest_column = strrep(CellTable.nearest_column,'_Barreloid','');
CellTable.layer = nan(size(CellTable.cortical_depth)); 

%% Add layer location of soma based on interpolated surfaces from
% depth values provided by Meyer et al.
tblLayer = readtable([inputPath 'neurons_layer.csv']);

for lay = [1:6]
   idxtmp = (tblLayer.layer==lay); 
   idstmp = tblLayer.id(idxtmp);
   idxtmp2 = ismember(CellTable.neuronID,idstmp);
   CellTable.layer(idxtmp2) = lay; 
end

%%
save([outputPath 'CellTable_' modelID '.mat'],'CellTable');