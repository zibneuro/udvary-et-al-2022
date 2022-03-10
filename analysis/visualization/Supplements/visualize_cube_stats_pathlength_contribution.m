%% Figure S2F
% Plot as boxplot for subvolume of 100x100x50 in C2 Column
% (to compare to Motta et al., Science 2019)
% only intracortical!
% - Axon contributing cells
% - Dendrite contributing cells
% - Axon contributing cell types
% - Dendrite contributing cell types
% - Axon path length to soma
% - Dendrite path length to soma
%
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
figurePath = [matlabPath 'output\Figures\']; 
tablePath = [matlabPath 'output\Tables\']; 
filepath = [matlabPath 'data\dataSubvolume\'];

%%  
f1 = figure(1);
clf;

% Axon Contributing Cells
load([filepath 'pre_CellTypeContribution_100_100_50.mat']);
axonCells = sum(ctContribution,2);
subplot(1,6,1);
boxplot(axonCells./1e3,'Symbol','k+','Color','k');
set(gca,'Box','off','TickDir','out','XTick',[], ...
    'Ytick',[0:5:50]);
title('axon #cells (1e3)');

% Dendrite Contributing Cells
load([filepath 'post_CellTypeContribution_100_100_50.mat']);
dendCells = sum(ctContribution,2);
subplot(1,6,2);
boxplot(dendCells./1e3,'Symbol','k+','Color','k');
set(gca,'Box','off','TickDir','out','XTick',[], ...
    'Ytick',[1.6:0.4:4]);
title('dend #cells (1e3)');

% Axon Contributing Cell Types
load([filepath 'pre_CellTypeContribution_100_100_50.mat']);
axonCellTypes = sum(ctContribution>0,2);
subplot(1,6,3);
boxplot(axonCellTypes,'Symbol','k+','Color','k');
set(gca,'Box','off','TickDir','out','XTick',[], ...
    'Ytick',[0:1:20]);
title('axon #celltypes');

% Dendrite Contributing Cell Types
load([filepath 'post_CellTypeContribution_100_100_50.mat']);
dendCellTypes = sum(ctContribution>0,2);
subplot(1,6,4);
boxplot(dendCellTypes,'Symbol','k+','Color','k');
set(gca,'Box','off','TickDir','out','XTick',[], ...
    'Ytick',[0:2:20]);
title('dend #celltypes');

% Axon Path Length to soma
load([filepath 'pre_pathLength_100_100_50.mat']);
axonPathLen = len;
subplot(1,6,5);
boxplot(axonPathLen./1e3,'Symbol','','Color','k');
hold on;
plot(1,max(axonPathLen./1e3),'k+');
set(gca,'Box','off','TickDir','out','XTick',[], ...
    'Ytick',[0:2:10]);
title('axon pathlen (mm)');

% Dendrite Path Length to soma
load([filepath 'post_pathLength_100_100_50.mat']);
dendPathLen = len;
subplot(1,6,6);
boxplot(dendPathLen./1e3,'Symbol','','Color','k');
hold on;
plot(1,max(dendPathLen./1e3),'k+');
set(gca,'Box','off','TickDir','out','XTick',[], ...
    'Ytick',[0:0.5:2]);
title('dend pathlen (mm)');

%% Save figure
set(f1,'PaperPositionMode','auto','Position',[0 0 950 200]); 
print(f1,'-painters','-dsvg','-r600',[figurePath 'CubeStats_2.svg']); 

%% Write Results in .csv file
optStringParam = {'axonPathLengthToSoma[mm]','dendPathLengthToSoma[mm]', ...
            'totalPathLengthToSoma[mm]','#axonCells','#dendCells', ...
            '#totalCells','#axonCellTypes','#dendCellTypes'};

fid = fopen([tablePath 'CubeStats_2_100x100x50.csv'],'w');
if fid==-1
    error('Cannot write file');
end
fprintf(fid,'Param,Mean,SD,Median,25th,75th,min,max,n,\n');

% - Axon contributing cells
% - Dendrite contributing cells
% - Axon contributing cell types
% - Dendrite contributing cell types
% - Axon path length to soma
% - Dendrite path length to soma
for optParam = 1:length(optStringParam)

    switch optParam
        case 1
           x = axonPathLen./1e3;
        case 2
           x = dendPathLen./1e3;
        case 3
           x = [dendPathLen; axonPathLen]./1e3;
        case 4
           x = axonCells;
        case 5
           x = dendCells; 
        case 6
           x = axonCells + dendCells;
        case 7
           x = axonCellTypes;
        case 8
           x = dendCellTypes;
    end

    fprintf(fid,'%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,\n', ...
                optStringParam{optParam},mean(x),std(x), ...
                median(x),prctile(x,25),prctile(x,75), ...
                min(x),max(x),numel(x));
end
fclose(fid);