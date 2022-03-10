%% Figure 3D
% example connection probability distribution
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
figurePath = [matlabPath 'output\Figures\']; 
tablePath = [matlabPath 'output\Tables\']; 
dataPath = [matlabPath 'data\example_selected_cells\'];
binsz = 0.02; 
XSample = [0:binsz:1+binsz];
            
%%
L5PT_ID = 301854;
L2PY_ID = 748854;
L6ACC_ID = 199678;

ontoExample = [];
fromExample = [];

for idtmp = [L5PT_ID L2PY_ID L6ACC_ID]
    preTbl = readtable([dataPath num2str(idtmp) '\' num2str(idtmp) ...
                        '_DSC.csv']);
    postTbl = readtable([dataPath num2str(idtmp) '\' num2str(idtmp) ...
                        '_innervating-pre-neurons-dsc.csv']);    
               
    pre_ids = idtmp.*ones(size(preTbl.DSC)); 
    fromExample = [fromExample; pre_ids preTbl.post_id preTbl.DSC];
                    
    post_ids = idtmp.*ones(size(postTbl.dsc)); 
    ontoExample = [ontoExample; postTbl.pre_id post_ids postTbl.dsc];
end

%% Generate matrix view
pre_ids = unique(ontoExample(:,1));
post_ids = unique(ontoExample(:,2)); 
numZeros = numel(pre_ids)*numel(post_ids)-size(ontoExample,1);
val = [ontoExample(:,3); zeros(numZeros,1)];

pre_ids = unique(fromExample(:,1));
post_ids = unique(fromExample(:,2)); 
numZeros = numel(pre_ids)*numel(post_ids)-size(fromExample,1);
val = [val; fromExample(:,3); zeros(numZeros,1)];
val = 1-exp(-val); 
valNoZero = val(val>0); 

%% With zeros
f1 = figure(1);
clf;

% Load Innervation values
N = histcounts(val,XSample);
maxVal = max(N); 

% Plot
stairs([XSample(1) XSample(1:end-1) XSample(end)].*100,[0 N 0], ...
            'Color','k'); 
hold on;
errorbar(mean(val).*100,1.08*maxVal,std(val).*100,'horizontal','ko', ...
    'MarkerEdgeColor','none','MarkerFaceColor','k');
ylabel('occurrences');
xlabel(['P(A,B) (%)']);
set(gca,'TickDir','out','Box','off','XLim',[-2 100],'XTick',[0:25:100], ...
    'YTick',[1 maxVal/2 maxVal],'YTickLabel','','Yscale','linear');  

set(f1,'PaperPositionMode','auto','Position',[0 0 250 550]); 
print(f1,'-dsvg','-r600',[figurePath 'ExampleCPDistr_Triplet.svg']); 

%% With zeros
f2 = figure(2);
clf;

% Load Innervation values
N = histcounts(valNoZero,XSample);
maxVal = max(N); 

% Plot
stairs([XSample(1) XSample(1:end-1) XSample(end)].*100,[0 N 0], ...
            'Color','k'); 
hold on;
errorbar(mean(valNoZero).*100,1.08*maxVal,std(valNoZero).*100,'horizontal','ko', ...
    'MarkerEdgeColor','none','MarkerFaceColor','k');
ylabel('occurrences');
xlabel(['P(A,B) (%)']);
set(gca,'TickDir','out','Box','off','XLim',[-2 100],'XTick',[0:25:100], ...
    'YTick',[1 maxVal/2 maxVal],'YTickLabel','','Yscale','linear');  

set(f2,'PaperPositionMode','auto','Position',[0 0 250 150]); 
print(f2,'-dsvg','-r600',[figurePath 'ExampleCPDistr_Triplet_NoZeros.svg']); 

tbl = table();
tbl.min_connection_probability = XSample(1:end-1)';
tbl.occurrences = N'; 
writetable(tbl,[tablePath 'ConnectionProbabilityHistogram_ExampleIllustration.csv']);