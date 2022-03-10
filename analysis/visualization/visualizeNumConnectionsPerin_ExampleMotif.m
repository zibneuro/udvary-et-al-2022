%% Figure 6F
% uses slice data
% Compare connected neurons per L5PT motifs (different motif sizes)
% vs. motif occurrences vs. random motif occurrences
% Comparison to Perin et al.,
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
figurePath = [matlabPath 'output\Figures\']; 
tablePath = [matlabPath 'output\Tables\']; 
storagepath = [matlabPath 'data\higherOrder\dataPerin\']; 

%%
numTrials = 1000;          % number of neurons extracted
samplingThreshold = 1000;   % number of random samples for each edge 

f1 = figure(1);
clf;
numNodes = 8;
numDirEdges = numNodes * numNodes - numNodes; 
strNumNodes = num2str(numNodes); 
sliceIDs = [0:1:19];

p1 = nan(length(sliceIDs),57);
p2 = nan(length(sliceIDs),57);

for i = 1:length(sliceIDs)
    load([storagepath strNumNodes '_' num2str(numTrials) '_Sampling_' ...
                    num2str(samplingThreshold) '_Slice-' ...
                    num2str(sliceIDs(i)) '.mat']); 
                
    % Correct for sampling: incorrect absolute probability values
    if sum(numCombinationsList>samplingThreshold)>0

        numCombinationsListAll = [numCombinationsList(1:end-1) ...
                                            fliplr(numCombinationsList)];
        idx = numCombinationsListAll>samplingThreshold; 

        tmp_pMotif_avg = pMotif_avg(idx).*numCombinationsListAll(idx);
        tmp_pMotifUniform = pMotifUniform(idx).*numCombinationsListAll(idx);

        corrFactor_pMotif = (1-sum(pMotif_avg(~idx)))/sum(tmp_pMotif_avg);
        corrFactor_pMotifUniform = (1-sum(pMotifUniform(~idx)))/sum(tmp_pMotifUniform);

        pMotif_avg(idx) = tmp_pMotif_avg.*corrFactor_pMotif;
        pMotifUniform(idx) = tmp_pMotifUniform.*corrFactor_pMotifUniform;

        fprintf('[%d]: CorrectionFactor = [%.2f %.2f] [%.2f %.2f]\n', ...
                        numNodes,corrFactor_pMotif, ...
                        corrFactor_pMotifUniform,sum(pMotif_avg), ...
                        sum(pMotifUniform));           
    end       
    
    p1(i,:) = pMotif_avg;
    p2(i,:) = pMotifUniform; 
end

pDev = p1./p2;
pDev = mean(pDev); % mean here can distort values (>1 vs <1)

pMotif_avg = mean(p1);
pMotifUniform = mean(p2); 

fprintf('[%d]: %d of %d [%.2f%%] are underrepresented\n', ...
            numNodes,sum(pDev<1),numel(pDev),sum(pDev<1)/numel(pDev)*100);
pDevNorm = pDev; %./max(pDev); 
pDevNorm(pDev<1) = nan; 

%     pDevNorm = pDevNorm-min(pDevNorm); % [0 to Inf]; 
pDevNorm = pDevNorm./max(pDevNorm); % [0 to 1]; 

fprintf('     [%.2e %.2e]\n',min(pDevNorm),max(pDevNorm));
fprintf('     SUM: [%.2f %.2f]\n',sum(pMotif_avg),sum(pMotifUniform));

xValues = linspace(0,1,numDirEdges+1);

%%
figure(1);
clf;
plot(0:numDirEdges,pMotif_avg,'k.-'); 
hold on;
plot(0:numDirEdges,pMotifUniform,'.-','Color',[0.5 0.5 0.5]); 
set(gca,'Box','off','TickDir','out','YScale','log', ...
    'XLim',[0 numDirEdges],'Xtick',[0:10:numDirEdges]); 
    
set(f1,'PaperPositionMode', 'auto','Position',[0 0 300 150]);
print(f1,'-dsvg','-r600',[figurePath 'Perin_NumConnectionsExample.svg']);

%% Save as table
t = table([0:numDirEdges]',pMotif_avg',pMotifUniform', ...
                'VariableNames',{'NumEdges','pMotif(Model)','pMotif(Random)'});
writetable(t,[tablePath 'Perin_NumConnections_8Motif.csv']);