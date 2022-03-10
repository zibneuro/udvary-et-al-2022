%% Figure 6A & 6B
% Connection per branch pair in 100x100x50um
% Comparison to Motta et al.
% Histogram of unconnected branch pairs
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
figurePath = [matlabPath 'output\Figures\']; 
tablePath = [matlabPath 'output\Tables\']; 
dataPath = [matlabPath 'data\connections-per-pair\100-100-50\'];
volVoxel = 100*100*50;

%% Measured data
% Motta et al.
dataMottaNumContacts = [1 2 3 4];
numSyn = 388554;
dataMottaFrequency = [numSyn-sum([6045 474 83]) 6045 474 83]; % number is from biorxiv version!
% Fig3: 388,554 synapses detected 
% (153,171 synapses with more than 10 synapses in volume) 
% Other synapses numbers (6045 474 83) from biorxiv version!
% https://www.biorxiv.org/content/10.1101/460618v1.article-info 
% see Line 478 on page 17 
volMotta = 61.8*94.8*92.6;
        
%% Branch Pairs
tbl = readtable([dataPath 'summary_branch.csv']);
numContacts = [1:4];
numOccurrences = [tbl.k1 tbl.k2 tbl.k3 tbl.k4+tbl.k5+ ...
                        tbl.k6+tbl.k7+tbl.k8+tbl.k_8];
percentageUnconnected = tbl.k0./(sum(numOccurrences,2)+tbl.k0).*100;           
medOcc = median(numOccurrences);
minOcc = min(numOccurrences);
maxOcc = max(numOccurrences); 

% Calcualte to percentages (for text)
fprintf('Unconnected branch pairs: %.2f +/- %.2f\n', ...
        mean(percentageUnconnected),std(percentageUnconnected));
t2 = table2array(tbl(:,4:end));
fprintf('  #pairs: %.2e\n',median(sum(t2,2)));
for i = 1:length(numContacts)
   fprintf('  %d: %.2e\n',numContacts(i),medOcc(i)); 
end

tblResult = table; 
tblResult.ConnectionsPerBranchPair = num2cell(numContacts'); 
tblResult.ConnectionsPerBranchPair(end) = {['>=' num2str(numContacts(end))]};
tblResult.medianOccurrences = medOcc';
tblResult.minOccurrences = minOcc';
tblResult.maxOccurrences = maxOcc';
tblResult.numVolumes = size(numOccurrences,1).*ones(size(tblResult.maxOccurrences));
writetable(tblResult,[tablePath 'ConnectionsPerBranchPair_L4_100x100x50.csv']); 

tblResult = table; 
tblResult.Subvolume = [1:numel(percentageUnconnected)]'; 
tblResult.PercentageUnconnected = percentageUnconnected; 
writetable(tblResult,[tablePath 'UnconnectedBranchPairs_L4_100x100x50.csv']); 

%%
f1 = figure(1);
clf;
plot(numContacts,medOcc,'k-');
hold on;
plot(numContacts,minOcc,'k--');
plot(numContacts,maxOcc,'k--');
plot(dataMottaNumContacts,dataMottaFrequency,'o','MarkerFaceColor','g', ...
        'MarkerEdgeColor','none');
xlabel('connections per branch pair');
ylabel('occurrences per 500,000 um3');
set(gca,'Box','off','TickDir','out','YScale','log','Xtick',[1:4], ...
    'YTick',[1e0 1e2 1e4 1e6 1e8],'Xlim',[1 4],'YLim',[5 1e6]);
axis square; 

f2 = figure(2);
h = histogram(percentageUnconnected,[0:0.02:100],'FaceColor','none', ...
        'EdgeColor','k','FaceAlpha',1,'DisplayStyle','stairs');
set(gca,'Box','off','TickDir','out','XLim',[99.7 99.9], ...
    'Ylim',[0 max(h.Values)],'YTick',[0 max(h.Values)],'Xtick',[99:0.05:100]);

%% Cell Pairs
tbl = readtable([dataPath 'summary_cell.csv']);           
numContacts = [1:4];
numOccurrences = [tbl.k1 tbl.k2 tbl.k3 tbl.k4+tbl.k5+ ...
                        tbl.k6+tbl.k7+tbl.k8+tbl.k_8];
percentageUnconnected = tbl.k0./(sum(numOccurrences,2)+tbl.k0).*100;                             
medOcc = median(numOccurrences);
minOcc = min(numOccurrences);
maxOcc = max(numOccurrences); 

% Calcualte to percentages (for text)
fprintf('Unconnected cell pairs: %.2f +/- %.2f\n', ...
        mean(percentageUnconnected),std(percentageUnconnected));
t2 = table2array(tbl(:,4:end));
fprintf('  #pairs: %.2e\n',median(sum(t2,2)));
for i = 1:length(numContacts)
   fprintf('  %d: %.2e\n',numContacts(i),medOcc(i)); 
end

fprintf('------------\nSample size = %d\n',size(numOccurrences,1));

%% Save figures
set(f1,'PaperPositionMode','auto','Position',[0 0 250 250]); 
print(f1,'-painters','-dsvg','-r600',[figurePath ...
                    'connection_per_branchpair_100_100_50_V2.svg']); 
set(f2,'PaperPositionMode','auto','Position',[0 0 250 150]); 
print(f2,'-painters','-dsvg','-r600',[figurePath ...
        'UnconnectedBranches_Motta.svg']); 