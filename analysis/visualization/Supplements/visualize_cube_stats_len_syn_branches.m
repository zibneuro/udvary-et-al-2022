%% Figure S2D
% Plot as boxplot for subvolume of 100x100x50 in C2 Column
% (to compare to Motta et al., Science 2019)
% - Dendrite Length
% - Axon Length
% - Number of Dendrite Branches
% - Number of Axon Branches
% - Number of (Axon) Boutons
%
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
figurePath = [matlabPath 'output\Figures\']; 
tablePath = [matlabPath 'output\Tables\']; 
filepath = [matlabPath 'data\cube-stats_ref-volume\'];

%%      
f1 = figure(1);
clf;

% Axon / Pre  
filename = 'cube-stats_100-100-50_ref-volume_pre.csv';
tbl = readtable([filepath filename]);

% Axon Length
subplot(1,5,1);
axonLen = tbl.length;
boxplot(axonLen/1e6,'Symbol','k+','Color','k');
set(gca,'Box','off','TickDir','out','XTick',[]);
title('axon length (m)');

% Axon Branches
subplot(1,5,3);
axonBranches = tbl.branches;
boxplot(axonBranches/1e3,'Symbol','k+','Color','k');
set(gca,'Box','off','TickDir','out','XTick',[], ...
    'YTick',[20:10:80]);
title('axon branches (1e3)');

% Axon Boutons (=Synapses)
subplot(1,5,5);
syn = tbl.boutons;
boxplot(syn/1e5,'Symbol','k+','Color','k');
set(gca,'Box','off','TickDir','out','XTick',[], ...
    'YTick',[2:2:10]);
title('synapses (1e5)');

% Dendrite / Post
filename = 'cube-stats_100-100-50_ref-volume_post.csv';
tbl = readtable([filepath filename]);

% Dendrite Length
subplot(1,5,2);
dendLen = tbl.length;
boxplot(dendLen/1e6,'Symbol','k+','Color','k');
set(gca,'Box','off','TickDir','out','XTick',[], ...
    'YTick',[0.2:0.05:0.4]);
title('dend length (m)');

% Dendrite Branches
subplot(1,5,4);
dendBranches = tbl.branches; 
boxplot(dendBranches/1e3,'Symbol','k+','Color','k');
set(gca,'Box','off','TickDir','out','XTick',[], ...
    'YTick',[4:1:9]);
title('dend branches (1e3)');

%% Save figure
set(f1,'PaperPositionMode','auto','Position',[0 0 800 200]); 
print(f1,'-painters','-dsvg','-r600',[figurePath 'CubeStats_1.svg']); 

%% Write Results in .csv file
optStringParam = {'axonLength[m]','dendLength[m]','totalLength[m]', ...
                '#axonBranchlets','#dendBranchlets','totalBranchlets', ...
                '#boutons'};

fid = fopen([tablePath 'CubeStats_1_100x100x50.csv'],'w');
if fid==-1
    error('Cannot write file');
end
fprintf(fid,'Param,Mean,SD,Median,25th,75th,min,max,n,\n');

% Axon/Dend Length
% Axon/Dend Branches
% AxonBoutons
for optParam = 1:length(optStringParam)

    switch optParam
        case 1
           x = axonLen./1e6;
        case 2
           x = dendLen./1e6;
        case 3
            x = dendLen./1e6 + axonLen./1e6;
        case 4
           x = axonBranches;
        case 5
           x = dendBranches; 
        case 6
            x = dendBranches + axonBranches;
        case 7
           x = syn;
    end

    fprintf(fid,'%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,\n', ...
                optStringParam{optParam},mean(x),std(x), ...
                median(x),prctile(x,25),prctile(x,75), ...
                min(x),max(x),numel(x));
end
fclose(fid);