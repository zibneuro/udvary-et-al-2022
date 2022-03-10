%% Figure 2E
% Connections per branch pair vs. occurrencees
% Connections per cell pair vs. occurrencees
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
figurePath = [matlabPath 'output\Figures\']; 
tablePath = [matlabPath 'output\Tables\']; 
dataPath = [matlabPath 'data\connections-per-pair\'];

%% Branch Pairs
tbl = readtable([dataPath 'summary_branch.csv']);
numContacts = [0:8];
numOccurrences = [tbl.k0 tbl.k1 tbl.k2 tbl.k3 tbl.k4 tbl.k5 ...
                        tbl.k6 tbl.k7 tbl.k8+tbl.k_8];
                    
medOcc = median(numOccurrences);
minOcc = min(numOccurrences);
maxOcc = max(numOccurrences); 

% Calcualte to percentages (for text)
percentageUnconnected = numOccurrences(:,1)./sum(numOccurrences,2).*100;
fprintf('Unconnected branch pairs (in %%): %.2f +/- %.2f\n', ...
        mean(percentageUnconnected),std(percentageUnconnected));
    
for i = 1:5
    percentageConnected = numOccurrences(:,numContacts==i)./ ...
                            sum(numOccurrences(:,numContacts>0),2).*100;
    fprintf('Connected branch pairs with %d synapses (in %%): %.4f +/- %.4f\n', ...
            i,mean(percentageConnected),std(percentageConnected));
end

tblResult = table; 
tblResult.ConnectionsPerBranchPair = num2cell(numContacts'); 
tblResult.ConnectionsPerBranchPair(end) = {['>=' num2str(numContacts(end))]};
tblResult.medianOccurrences = medOcc';
tblResult.minOccurrences = minOcc';
tblResult.maxOccurrences = maxOcc';
tblResult.numVolumes = size(numOccurrences,1).*ones(size(tblResult.maxOccurrences));
writetable(tblResult,[tablePath 'ConnectionsPerBranchPair.csv']); 

%%
f1 = figure(1);
clf;
plot(numContacts,medOcc,'k-');
hold on;
plot(numContacts,minOcc,'k--');
plot(numContacts,maxOcc,'k--');
xlabel('connections per branch pair');
ylabel('occurrences per (50um)3');
set(gca,'Box','off','TickDir','out','YScale','log','Xtick',[0:2:8], ...
    'YTick',[1e0 1e2 1e4 1e6 1e8],'Xlim',[0 8],'YLim',[0.5 1e8]);

%% Cell Pairs
tbl = readtable([dataPath 'summary_cell.csv']);
numContacts = [0:8];
numOccurrences = [tbl.k0 tbl.k1 tbl.k2 tbl.k3 tbl.k4 tbl.k5 ...
                        tbl.k6 tbl.k7 tbl.k8+tbl.k_8];
                    
medOcc = median(numOccurrences);
minOcc = min(numOccurrences);
maxOcc = max(numOccurrences); 

% Calcualte to percentages (for text)
percentageUnconnected = numOccurrences(:,1)./sum(numOccurrences,2).*100;
fprintf('Unconnected cell pairs (in %%): %.2f +/- %.2f\n', ...
        mean(percentageUnconnected),std(percentageUnconnected));
    
percentageConnected5 = numOccurrences(:,numContacts==5)./ ...
                        sum(numOccurrences(:,numContacts>0),2).*100;
fprintf('Connected cell pairs with 5 synapses (in %%): %.2f +/- %.2f\n', ...
        mean(percentageConnected5),std(percentageConnected5));

%%
f2 = figure(2);
clf;
plot(numContacts,medOcc,'k-');
hold on;
plot(numContacts,minOcc,'k--');
plot(numContacts,maxOcc,'k--');
xlabel('connections per cell pair');
ylabel('occurrences per (50um)3');
set(gca,'Box','off','TickDir','out','YScale','log','Xtick',[0:2:8], ...
    'YTick',[1e0 1e2 1e4 1e6 1e8],'Xlim',[0 8],'YLim',[0.5 1e8]);

%% Save figures
set(f1,'PaperPositionMode','auto','Position',[0 0 400 250]); 
print(f1,'-painters','-dsvg','-r600',[figurePath ...
                    'connection_per_branchpair.svg']); 
set(f2,'PaperPositionMode','auto','Position',[0 0 400 250]); 
print(f2,'-painters','-dsvg','-r600',[figurePath ...
                    'connection_per_cellpair.svg']); 
saveas(f1,[figurePath 'connection_per_branchpair.jpg']);
saveas(f2,[figurePath 'connection_per_cellpair.jpg']);