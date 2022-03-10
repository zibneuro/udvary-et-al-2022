%% Generates 
% - AxonLenData.mat
% - DendriteLenDataAligned.mat
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
dataFolder = [matlabPath 'preprocessing\data\robustness_morphologies\'];

%% Axons
CellTypes = generateCellTypeList('excitatory'); 
CellTypes{end+1} = 'VPM'; 
CellTypeFolders = {'L2','L3-4','L4Py','L4SP','L4SS','L5ST','L5TT', ...
                        'L6CC','L6CC_inv','L6CT','VPM'};

lenMatrix = cell(1,length(CellTypes)); 
numSamples = nan(1,length(CellTypes)); 

for ct = 1:length(CellTypes)
    
    d = dir([dataFolder 'D2_Morphologies\Axons\' CellTypeFolders{ct} ...
            'axon\*_Axon_Density.am']);    
          
    numSamples(ct) = length(d);
    for i = 1:length(d)
        filename = [d(i).name(1:end-3) '.mat']; 
        filename = strrep(filename,'Axon_Density','AxonDensity'); 
        
        % Loads variables: BBx BBy BBz cellType matrix
        m = load([dataFolder 'dataDensity\' filename],'BBx','BBy','BBz','matrix');
        
        if i==1
            matrix = nan([length(d) size(m.matrix)]);
        end
        matrix(i,:,:,:) = m.matrix; 
    end
    lenMatrix{ct} = matrix; 
end

BBx = m.BBx;
BBy = m.BBy;
BBz = m.BBz;

save([dataFolder 'AxonLenData.mat'], ...
    'BBx','BBy','BBz','lenMatrix','CellTypes','numSamples');  

%% Dendrites Soma Aligned
CellTypes = generateCellTypeList('excitatory'); 
CellTypeFolders = {'L2','L3-4','L4Py','L4SP','L4SS','L5ST','L5TT', ...
                        'L6CC','L6CC_inv','L6CT'};
lenMatrix = cell(1,length(CellTypes)); 
numSamples = nan(1,length(CellTypes)); 

for ct = 1:length(CellTypes)
    
    d = dir([dataFolder 'D2_Morphologies\Dendrites\' CellTypeFolders{ct} ...
            '\*_BasalDendrite_Density.am']);    
          
    numSamples(ct) = length(d);
    for i = 1:length(d)
        filename = [d(i).name(1:end-3) '.mat']; 
        filename = strrep(filename,'BasalDendrite_Density','DendriteDensity'); 
        
        % Loads variables: BBx BBy BBz cellType matrix
        m = load([dataFolder 'dataDensityAligned\' filename],'BBx','BBy','BBz','matrix');
        
        if i==1
            matrix = nan([length(d) size(m.matrix)]);
        end
        matrix(i,:,:,:) = m.matrix; 
    end
    lenMatrix{ct} = matrix; 
end

BBx = m.BBx;
BBy = m.BBy;
BBz = m.BBz;

save([dataFolder 'DendriteLenDataAligned.mat'], ...
    'BBx','BBy','BBz','lenMatrix','CellTypes','numSamples'); 