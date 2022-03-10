function [matrixColor] = getBlueRedColorMatrix(matrix,switchValue)
%% Function that returns a 3D Matrix with color values for 2D Matrix
% in colormap [blue (=min) -> white (=switchValue) -> red (=max)]
%   Seperate linear spacing below and above switchValue (128 values each)
%   Both spacings might differ! Not symmetric unless 
%   -min(matrix(:)) == max(matrix(:))
% Input:
% - matrix: [n x m] matrix with values postivie and negative values
% - switchValue (optional; default: 0): value of the color white, 
%       values below will be bluer, values above will be reder
% Output:
% - matrixColor: [n x m x 3] matrix with color values
%       -> imagesc(1:size(matrixColor,1),1:size(matrixColor,2),matrixColor)
%
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
% Date: Feb, 9th, 2022

    if length(size(matrix))~=2
       error('Only works for 2D matrices');  
    end

    if nargin<2
        switchValue = 0;
    end
    
    colRed = [];
    xRed = [];
    colBlue = [];
    xBlue = [];
    
    if max(matrix(:))>switchValue    
        % Red Color Map for positive Values
        colRed = ones(128,3);
        xRed = linspace(switchValue,max(matrix(:)),size(colRed,1)); 
        colRed(:,2) = 1-linspace(0,1,size(colRed,1));
        colRed(:,3) = 1-linspace(0,1,size(colRed,1));
    end

    if min(matrix(:))<switchValue
        % Blue Color Map for negative Values
        colBlue = ones(128,3);
        xBlue = linspace(min(matrix(:)),switchValue,size(colBlue,1)); 
        colBlue(:,1) = linspace(0,1,size(colBlue,1));
        colBlue(:,2) = linspace(0,1,size(colBlue,1));
    end
    
    % Merge Colormaps
    col = [colBlue; colRed];
    x = [xBlue xRed];

    % 3D Matrix with Color Values
    matrixColor = nan(size(matrix,1),size(matrix,2),3);

    for i = 1:size(matrixColor,1) 
        for j = 1:size(matrixColor,2)
            binID = find(histcounts(matrix(i,j),x)==1);
            
            if isempty(binID)
                continue;
            end
            
            matrixColor(i,j,:) = col(binID,:);
        end
    end

end