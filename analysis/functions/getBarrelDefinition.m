function [BarrelCenter,BarrelzAxis,BarrelRadius] = ...
                    getBarrelDefinition(BarrelID)
% [BarrelCenter,BarrelzAxis,BarrelRadius] = getBarrelDefinition(BarrelID)
% Returns Barrel Definition like 3D center location of barrel, 
% its z-axis unit vector and its barrel radius
% Input:
% - filename: BarrelID (one of (Alpha, A1, A2, ... E3, E4))
% Output:
% - BarrelCenter: 3D center location of barrel [1x3]
% - BarrelzAxis: z-axis unit vector of barrel [1x3]
% - BarrelRadius: radius of barrel [1x1]
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)

    [BarrelIDs,BarrelCenter,BarrelzAxis,BarrelRadius] = ...
                                    readBarrelDefinitionSpreadSheet();
                                
    idx = strcmp(BarrelIDs,BarrelID);
    if sum(idx)~=1
       error('BarrelID %s does not exist!',BarrelID); 
    end
    
    BarrelCenter(~idx,:) = [];
    BarrelzAxis(~idx,:) = [];
    BarrelRadius(~idx) = [];
end