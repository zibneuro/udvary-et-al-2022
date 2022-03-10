%% Display mean displacement of axons during model generation
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
load([matlabPath 'data\CellTable_RBC_20.mat'],'CellTable');
tbl = readtable([matlabPath 'data\axon_displacement.csv']); 

dist = sqrt((tbl.x0-tbl.x1).^2+(tbl.y0-tbl.y1).^2+(tbl.z0-tbl.z1).^2); 
idxC2 = strcmp(CellTable.nearest_column,'C2');
idxEXC = ~strcmp(CellTable.cell_type,'INH') & ...
            ~strcmp(CellTable.cell_type,'VPM');
neuronID_vS1 = CellTable.neuronID(idxEXC); 
neuronID_C2 = CellTable.neuronID(idxC2 & idxEXC); 
clear CellTable;
idxC2 = ismember(tbl.neuron_id,neuronID_C2);
idxvS1 = ismember(tbl.neuron_id,neuronID_vS1);

%%
figure(1);
clf;
subplot(2,1,1);
histogram(dist(idxvS1));
set(gca,'Box','off','TickDir','out');
xlabel('displacement within vS1 (\mum)');
ylabel('frequency');
subplot(2,1,2);
histogram(dist(idxC2));
set(gca,'Box','off','TickDir','out');
xlabel('displacement within C2 (\mum)');
ylabel('frequency');

fprintf('Displacement(vS1): %.0f +/- %.0f um\n',mean(dist(idxvS1)), ...
                std(dist(idxvS1)));
fprintf('Displacement(C2): %.0f +/- %.0f um\n',mean(dist(idxC2)), ...
                std(dist(idxC2)));