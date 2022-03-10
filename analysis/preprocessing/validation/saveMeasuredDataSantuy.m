%% Extract synapse densities measured in samples across all layers from
% P14 rat hindlimb published in Santuy et al., Brain Structure and Function 
% 2018; Supplementary Table 1
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
outputPath = [matlabPath 'data\']; 
inputPath = [matlabPath 'preprocessing\data\']; 

[~, ~, raw] = xlsread([inputPath 'Measured_SantuyEtAl_BrainStructFunc_2018.xlsx'], ...
                    'Sheet1','A2:E30');

data = reshape([raw{:}],size(raw));
Layer = data(:,1);
ExcSyn = data(:,4);
InhSyn = data(:,5);

save([outputPath 'Measured_SantuyEtAl.mat'],'Layer','ExcSyn','InhSyn');