% [] = writeSurfacesBasedOnDensity(matrix,outputPrefix,
%                                   BoundingBox,VoxelSize,threshold)
% Writes ImageDataVolume (matrix) as surfaces
% Input:
% - matrix: 3D matrix 
% - outputPrefix: path/to/outputfilePrefix
% - BoundingBox: [1 x 3] vector with 3 elements of origin in indices 
% - VoxelSize: [1 x 1] Size of voxel
%       (default: 50) 
% - threshold: [1 x 2] Threshold for writing out surface
%       (default: [0 Inf]) 
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
function [] = writeSurfacesBasedOnDensity(matrix,outputPrefix, ...
                        BoundingBox,VoxelSize,threshold)

    if nargin<5
        threshold = [0 inf];
    end
    if nargin<4
       VoxelSize = 50; 
    end
                    
    BBx = BoundingBox(1):VoxelSize:BoundingBox(2);
    BBy = BoundingBox(3):VoxelSize:BoundingBox(4);
    BBz = BoundingBox(5):VoxelSize:BoundingBox(6);

    for x = 1:length(BBx)
        for y = 1:length(BBy)
            for z = 1:length(BBz)
                
                if matrix(x,y,z)>threshold(1) && matrix(x,y,z)<threshold(2)
                    BBmin = [BBx(x) BBy(y) BBz(z)]-VoxelSize/2;
                    BBmax = [BBx(x) BBy(y) BBz(z)]+VoxelSize/2;

                    str = sprintf('%d_%d_%d',BBx(x),BBy(y),BBz(z));

                    outputfilename = [outputPrefix str];
                    writeBoxAmiraSurface(BBmin,BBmax,outputfilename); 
                end
            end
        end
    end
end