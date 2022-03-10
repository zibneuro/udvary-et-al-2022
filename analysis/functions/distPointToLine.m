%% Distance of point to line
% Input:
% - pt: 2D or 3D point vector
% - v1: 2D or 3D point vector of line
% - v2: 2D or 3D point vector of line
% v1 and v2 define the line
% Output:
% - dist: distance of point pt to line vector v1 and v2
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
function dist = distPointToLine(pt, v1, v2)

    % If 2D, make it 3D by adding 0
    if (length(v1)==2 && length(v2)==2 && length(pt)==2)
       v1 = [v1 0];  
       v2 = [v2 0];  
       pt = [pt 0];      
    end
    
    if (length(v1)~=3 || length(v2)~=3 || length(pt)~=3)
        error('Vectors need to be same size. Only works for size 2 or 3!');
    end

    % all need to be row vectors!
    if iscolumn(v1)
        v1 = v1';
    end
    if iscolumn(v2)
        v2 = v2';
    end
    if iscolumn(pt)
        pt = pt';
    end   
    
    a = v1 - v2;
    b = pt - v2;
    dist = norm(cross(a,b)) / norm(a);
end