%% Merge Amira Surfaces
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
function [pt] = mergeSurface(filenames,outputfilename,cube_size)

    if nargin<3
        cube_size = 50;
    end

    pt = cell(0);
    Vertices = [];
    Triangles = [];
    for i = 1:length(filenames)
        
        [Vertices1,Triangles1] = readAmiraSurface(filenames{i});
        vertexIDOffset = size(Vertices,2);
        Vertices = [Vertices Vertices1];
        Triangles = [Triangles Triangles1+vertexIDOffset];
        
        % Get lines around each box
        BBmin = min(Vertices1,[],2); 
        x = BBmin(1) + [0 cube_size];
        y = BBmin(2) + [0 cube_size];
        z = BBmin(3) + [0 cube_size];
        pt{end+1} = [x(1) y(1) z(1); x(2) y(1) z(1); ...
                        x(2) y(2) z(1); x(1) y(2) z(1); ...
                        x(1) y(2) z(2); x(2) y(2) z(2); ...
                        x(2) y(1) z(2); x(1) y(1) z(2); ...
                        x(1) y(1) z(1)];  
        pt{end+1} = [x(1) y(1) z(1); x(1) y(2) z(1)];          
        pt{end+1} = [x(1) y(1) z(2); x(1) y(2) z(2)];
        pt{end+1} = [x(2) y(2) z(1); x(2) y(2) z(2)];   
        pt{end+1} = [x(2) y(1) z(1); x(2) y(1) z(2)];  
    end  
    % Deletes identical surface triangle
    [v,~,ic] = unique(Vertices','rows');
    t = ic(Triangles);
    writeAmiraSurface(v',t,outputfilename);
end