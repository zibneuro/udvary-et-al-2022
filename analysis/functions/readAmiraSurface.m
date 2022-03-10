function [Vertices,Triangles] = readAmiraSurface(filename)
% [Vertices,Triangles] = readAmiraSurface(filename)
% Reads Amira Surface
% Input:
% - fiilename: path/to/filename.surf
% Output:
% - Vertices: [3 x NumVertices] 3D coordinates of Vertices
% - Triangles: [3 x NumTriangles] ID of Vertices
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
    fileID = fopen(filename,'r');

    if fileID==-1
       error(['Cannot read AmiraSurfacefile: ' filename]); 
    end
    
    tline = fgetl(fileID); 
    
    foundVertices = 0;
    foundTriangles = 0; 
    
    while ischar(tline)
        
        idx = strfind(tline,'Vertices'); 
        if ~isempty(idx) && foundVertices==0
            foundVertices=1;
            numVertices = cell2mat(textscan(tline,'%*s %d','delimiter',' '));
        end
        
        idx = strfind(tline,'Triangles'); 
        if ~isempty(idx) && foundTriangles==0
            foundTriangles=1;
            numTriangles = cell2mat(textscan(tline,'%*s %d','delimiter',' '));
        end
        
        if foundTriangles==1
            Triangles = nan(3,numTriangles); 
            
            for i = 1:numTriangles
                tline = fgetl(fileID);
                C = textscan(tline,'%d','delimiter','');
                Triangles(:,i) = cell2mat(C);
            end
            
            foundTriangles = 0;
        end
        
        if foundVertices==1
            Vertices = nan(3,numVertices);
            
            for i = 1:numVertices
                tline = fgetl(fileID);
                C = textscan(tline,'%f','delimiter','');
                Vertices(:,i) = cell2mat(C); 
            end
            
            foundVertices = 2;
        end
        
        tline = fgetl(fileID); 
    end
    
    fclose(fileID); 
end