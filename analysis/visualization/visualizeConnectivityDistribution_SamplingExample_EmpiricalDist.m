%% Figure 4D (bottom)
% Histogram of L5PT to L5PT for different 3D intersomatic Distances
% [0 100] [100 200] and [200 300]
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
figurePath = [matlabPath 'output\Figures\']; 
tablePath = [matlabPath 'output\Tables\']; 
binsz = 0.02; 
XSample = [0:binsz:1+binsz];
load([matlabPath 'data\CellMatrix_C2\L5PT_L5PT.mat']);
load([matlabPath 'data\CellTable_RBC_20.mat'],'CellTable');
idxtmp = strcmp(CellTable.cell_type,'L5PT') & ...
            strcmp(CellTable.nearest_column,'C2');
CellTable(~idxtmp,:) = [];
neuronID = CellTable.neuronID;
somaPos = [CellTable.soma_x CellTable.soma_y CellTable.soma_z]; 
clear CellTable;

%% Compute distance matrix
n = size(somaPos,1);
distanceMatrix = nan(n,n); 
for i = 1:n
    for j = 1:n
        if i~=j
            distanceMatrix(i,j) = sqrt(sum((somaPos(i,:) ...
                                - somaPos(j,:)).^2));
        end
    end
end
    
%%
f1 = figure(1);
clf;
maxVal = 0; 
                    
for ii = 1:3
    
    switch ii
        case 1
            dist = [0 100];
            col = [0 1 0];
        case 2
            dist = [100 200]; 
            col = [1 0.5 0];
        case 3
            dist = [200 300];
            col = [0 0 1];
    end
    
    idxtmp = distanceMatrix>=dist(1) & distanceMatrix<=dist(2);
    
    % Load Innervation values
    p = 1-exp(-I.I(idxtmp));
    N = histcounts(p(:),XSample);
    maxVal = max([maxVal N]); 
    
    % Plot
    hold on;
    stairs([XSample(1) XSample(1:end-1) XSample(end)].*100,[0 N 0], ...
                'Color',col); 
    ylabel('occurrences');
    xlabel(['P(A,B) (%)']);
    set(gca,'TickDir','out','Box','off','XLim',[-2 100],'XTick',[0:25:100], ...
        'YLim',[0 maxVal],'YTick',[0 maxVal/2 maxVal],'YTickLabel','');  
    
    tbl = table();
    tbl.connection_probability = round(p,4); 
    writetable(tbl,[tablePath 'L5PT_L5PT_' num2str(dist(1)) '_' ...
                    num2str(dist(2)) '_connection_probabilities.csv']);
end

%% Save figure
set(f1,'PaperPositionMode','auto','Position',[0 0 250 150]); 
print(f1,'-dsvg','-r600',[figurePath 'ExampleDistributions_Dist_EmpDist.svg']); 