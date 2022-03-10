% [matrix, origin, VoxelSize, BoundingBox] = readAmira(filename)
% Reads Simple 3D ImageDataVolume/ScalarField Amira file.
% Returns matrix, coorindates of origin, VoxelSize, and BoundingBox
% Input:
% - filename: path/to/filename.am
% Output:
% - matrix: [N1 x N2 x N3] 3D matrix containing values of ImageDataVolume
% - origin: [1 x 3] relativ position of the origin in N1 N2 N3
% - VoxelSize: [1 x 3] size of each voxel (typically 50 x 50 x 50)
% - BoundingBox: [1 x 6] BoundingBox of ImageDataVolume. 
%       [xmin xmax ymin ymax zmin zmax] 3D centers of respective voxels
%       x_pos = [BoundingBox(1):VoxelSize(1):BoundingBox(2)]
%       y_pos = [BoundingBox(3):VoxelSize(2):BoundingBox(4)]
%       z_pos = [BoundingBox(5):VoxelSize(3):BoundingBox(6)]
%           -> center of each voxel
%       actual min and max position is x_pos -/+ VoxelSize(1)/2
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
function [matrix, origin, VoxelSize, BoundingBox] = readAmira(filename)

    % Open the text file.
    if (~(strcmp(filename(end-2:end), '.am')))
        filename = [filename '.am'];
    end
    
    fileID = fopen(filename,'r');
    
    if fileID==-1
       error(['Cannot read Amira File ' filename]); 
    end
    
    % Compute dimensionality and BoundingBox of input matrix
    dim = nan(1,3); 
    BoundingBox = nan(1,6); 
    tline = fgetl(fileID); 
    while ischar(tline) && (sum(isnan(dim)>0) || sum(isnan(BoundingBox)>0))
        
        idx = strfind(tline,'define Lattice'); 
        if ~isempty(idx)
            C = textscan(tline(idx:end),'%*s %*s %d %d %d','delimiter',' ');
            dim = cell2mat(C);
        end
        
        idx = strfind(tline,'BoundingBox'); 
        if ~isempty(idx)
            C = textscan(tline(idx:end),'%*s %d %d %d %d %d %d','delimiter',' ');
            BoundingBox = cell2mat(C);
        end
        tline = fgetl(fileID); 
    end
    
    if sum(isnan(dim)>0)
       error('Definition of Lattice could not be found!');  
    end
    
    if sum(isnan(BoundingBox)>0)
       error('Definition of BoundingBox could not be found!');  
    end
    
    BoundingBox = double(BoundingBox);     
    matrix = zeros(dim); 
    idx = 0; 
    
    % Fill matrix with the corresponding values
    while (ischar(tline) && idx <= numel(matrix))
        
        if strcmp(tline,'@1')
            idx = idx+1;    
            tline = fgetl(fileID); 
        end
        
        if (idx > 0 && idx <= numel(matrix))
            C = textscan(tline,'%f64');
            matrix(idx) = cell2mat(C); 
            idx = idx + 1; 
        end
        tline = fgetl(fileID); 
    end

    while strcmp(tline,'')
        tline = fgetl(fileID); 
    end
    
    % Close the text file.
    fclose(fileID);
    
    if ischar(tline)
       warning(['Something did not work! There are still elements in the list that should be inserted into the matrix. ' tline]); 
       
    end
    
    % Compute VoxelSize
    VoxelSize(1) = abs(BoundingBox(1)-BoundingBox(2))/(dim(1)-1); 
    VoxelSize(2) = abs(BoundingBox(3)-BoundingBox(4))/(dim(2)-1); 
    VoxelSize(3) = abs(BoundingBox(5)-BoundingBox(6))/(dim(3)-1); 
    
    if VoxelSize(1) ~= VoxelSize(2) || VoxelSize(2) ~= VoxelSize(3)
       warning(['VoxelSize is not the same for all dimensionalities! ' ...
           num2str(VoxelSize)]) 
    end
    VoxelSize = double(VoxelSize); 

    % Compute Origin
    origin = BoundingBox([1 3 5])./VoxelSize; 
    origin = abs(origin) + 1; 
end