function [BarrelID,BarrelCenter,BarrelzAxis,BarrelRadius] = ...
                    readBarrelDefinitionSpreadSheet(filename)
% [BarrelID,BarrelCenter,BarrelzAxis,BarrelRadius] = ...
%                     readBarrelDefinitionSpreadSheet(filename)
% Reads in BarrelDefinitionSpreadSheet (for NeuroNet) and extracts parameters
% like 3D center location of each barrel, their z-axis unit vector and
% their barrel radius
% Input:
% - filename: path/to/BarrelDefinitionSpreadsheet.csv (optional)
% Output:
% - BarrelID: Label of Barrel [24x1] (Alpha, A1, A2, ... E3, E4)
% - BarrelCenter: 3D center location of each barrel [24x3]
% - BarrelzAxis: z-axis unit vector of each barrel [24x3]
% - BarrelRadius: radius of each barrel [24x1]
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
    if nargin<1
       filename = ['D:\udvary-et-al-2022\analysis\' ...
                    'preprocessing\data\' ...
                    'BarrelDefinitionSpreadsheet.csv'];
    end

    % Read in file
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, '%s%s%s%f%[^\n\r]','Delimiter','\t', ...
        'HeaderLines',1,'ReturnOnError',false,'EndOfLine','\r\n');
    fclose(fileID);

    % Formatting (get rid of "")
    BarrelID = strrep(dataArray{:,1},'"','');
    cleanCenter = strrep(dataArray{:,2},'"','');
    cleanZAxis = strrep(dataArray{:,3},'"','');

    % Extract features
    BarrelRadius = dataArray{:,4};
    BarrelCenter = nan(length(BarrelRadius),3);
    BarrelzAxis = nan(length(BarrelRadius),3); 

    for i = 1:length(BarrelRadius)
        BarrelCenter(i,:) = cell2mat(textscan(cleanCenter{i},'%f%f%f','Delimiter',','));
        BarrelzAxis(i,:) = cell2mat(textscan(cleanZAxis{i},'%f%f%f','Delimiter',','));
    end
    
end