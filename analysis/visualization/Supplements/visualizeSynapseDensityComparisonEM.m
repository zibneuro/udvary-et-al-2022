%% Figure S2E
% plot Boxplot of synapse density in each layer of D2 barrel
% and compare to measured synapse density
% by Santuy et al., 2018 and Kasthuri et al., 2015
% Save 2 figures one for exc and inh synapses
% Save .csv file with results of comparison (only Santuy)
%
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
addpath([matlabPath 'functions\']);
load([matlabPath 'data\SynapseDensitiesAlongDepth.mat']);
figurePath = [matlabPath 'output\Figures\']; 
tablePath = [matlabPath 'output\Tables\']; 
colGreen = [177 200 0]./255; 
LayerBordersD2 = [0 157 296 575 900 1411 1973];
boolAll = 1; 

%% Measured synapse densities in different layers 
% Santuy, A., Rodriguez, J. R., DeFelipe, J., & Merchan-Perez, A. (2018). 
% Volume electron microscopy of the distribution of synapses in the 
% neuropil of the juvenile rat somatosensory cortex. 
% Brain Structure and Function, 223(1), 77?90.
% From Supplementary Table 1
measuredSantuy = load([matlabPath 'data\Measured_SantuyEtAl.mat']);  
tblMeasured = table();
tblMeasured.excitatory_synapse_density = measuredSantuy.ExcSyn;
tblMeasured.inhibitory_synapse_density = measuredSantuy.InhSyn;
tblMeasured.layer = measuredSantuy.Layer;
writetable(tblMeasured,[tablePath 'SynapseDensities_Measured.csv']);

%%
if boolAll
    voxelDensity = synDensity.all_density;
    voxelDepth = synDensity.all_depth;
else
    voxelDensity = synDensity.exc_density;
    voxelDepth = synDensity.exc_depth;
end

gModel = [];
synModel = []; 
gMeasured = [];
synMeasured = [];

% Go through all six layers
for i = 1:length(LayerBordersD2)-1

    % Measured values
    idxMeasured = (measuredSantuy.Layer==i);
    
    synMeasured_tmp = measuredSantuy.ExcSyn(idxMeasured)'; 
    
    if boolAll
        synMeasured_tmp = synMeasured_tmp + ...
                        measuredSantuy.InhSyn(idxMeasured)'; 
    end
    synMeasured = [synMeasured synMeasured_tmp];
    gMeasured = [gMeasured i.*ones(1,numel(synMeasured_tmp))]; 
    
    % Model's estimate
    idxtmp = LayerBordersD2(i)<=voxelDepth & ...
                LayerBordersD2(i+1)>=voxelDepth;
            
    if (sum(idxtmp)==0)
       error('No samples in layer found!'); 
    end
    
    synModel_tmp = voxelDensity(idxtmp);
    synModel = [synModel synModel_tmp];
    gModel = [gModel i.*ones(1,numel(synModel_tmp))]; 

    % Display Results
    fprintf('Layer %d:\n',i);
    synMeasuredDisp = synMeasured_tmp; 
    strTxt = ['Santuy et al., 2018'];
    fprintf(' Model: M = %.2f STD = %.2f Range = [%.2f %.2f] n = %d\n', ...
        mean(synModel_tmp),std(synModel_tmp), ...
        min(synModel_tmp),max(synModel_tmp),numel(synModel_tmp));
    fprintf(' Measured: M = %.2f STD = %.2f Range = [%.2f %.2f] n = %d (%s)\n', ...
        mean(synMeasuredDisp),std(synMeasuredDisp), ...
        min(synMeasuredDisp),max(synMeasuredDisp),numel(synMeasuredDisp), ...
        strTxt);
    
    % KS test
    [h,p,ks2stat] = kstest2(synModel_tmp,synMeasured_tmp,'Alpha',0.01);
    
    % h = 1 not from same distribution
    if h == 1
       fprintf(['synapses: Layer %d different distributions ' ...
           '[Santuy]: p = %.4f ks2stat = %.4f\n'],i,p,ks2stat); 
    end
end

%% Display Figure
f1 = figure(1);
clf;
boxplot(synModel,gModel,'Colors','k','Symbol','k+');
hold on;
plot(gMeasured+0.5,synMeasured,'o','MarkerEdgeColor','none', ...
            'MarkerFaceColor',colGreen,'MarkerSize',10,'LineWidth',1.5);     
set(gca,'Box','off','TickDir','out','XTick',[1:6]);
xlim([0.5 6.5]);
ylim([0 1.02*max(synModel)]);
xlabel('layer');
ylabel('synapse density per um');
set(f1,'PaperPositionMode', 'auto','Position',[0 0 350 250]);

if boolAll
    print(f1,'-dsvg','-r600',[figurePath 'AllSynapseDensity_BoxPlot.svg']); 
    t = table(synModel',gModel','VariableNames',{'AllBoutonDensity','Layer'});
    writetable(t,[tablePath 'EstimatedSynapseDensities_All.csv']);
else
    print(f1,'-dsvg','-r600',[figurePath 'ExcSynapseDensity_BoxPlot.svg']); 
    t = table(synModel',gModel','VariableNames',{'ExcBoutonDensity','Layer'});
    writetable(t,[tablePath 'EstimatedSynapseDensities_Exc.csv']);
end

%% Table with results
fprintf('Layer,Measured,PredictedMean,PredictedSD,PredictedMin,PredictedMax,\n');

for i = 1:6
   predictedSyn = synModel(i==gModel);  
   measuredSyn = synMeasured(i==gMeasured);  
   
   for j = 1:length(measuredSyn)
       fprintf('L%d,',i);
       fprintf('%.4f,',measuredSyn(j), ...
            mean(predictedSyn),std(predictedSyn),min(predictedSyn),max(predictedSyn));
       fprintf('\n');
   end
end