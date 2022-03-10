function [] = writeBoxAmiraSurface(BBmin,BBmax,outputfilename)
% [] = writeBoxAmiraSurface(BBmin,BBmax,outputfilename)
% Write Amira Surface 
% Input:
% - BBmin: [1 x 3] minimum coordinates of box
% - BBmax: [1 x 3] maximum coordinates of box
% - outputFilename: path/to/outputfilename.surf
% NOTE: have a look at the more flexible function writeAmiraSurface.m
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
    if (strcmp(outputfilename(end-4:end), '.surf'))
        fname = outputfilename;  
    else
        fname = [outputfilename '.surf'];
    end

    fid = fopen(fname,'w');
    
    if fid==-1
       error('Cannot write to %s',outputfilename); 
    end
    
    % Write Header
    fprintf(fid,['# HyperSurface 0.1 ASCII\n\n' ...
        'Parameters {\n' ...
        '    Materials {\n' ...
        '        Exterior {\n' ...
        '            id 0,\n' ...
        '            Color 1 1 1\n' ...
        '        }\n' ...
        '         {\n' ...
        '            id 1,\n' ...
        '            Color 0 0 0\n' ...
        '        }\n' ...
        '    }\n' ...
        '    BoundaryIds {\n' ...
        '        name "BoundaryConditions"\n' ...
        '    }\n' ...
        '}\n\n' ...
        'Vertices 8\n']);

    fprintf(fid,'\t%f %f %f\n',BBmax(1),BBmin(2),BBmin(3));
    fprintf(fid,'\t%f %f %f\n',BBmax(1),BBmax(2),BBmin(3));
    fprintf(fid,'\t%f %f %f\n',BBmin(1),BBmax(2),BBmin(3));
    fprintf(fid,'\t%f %f %f\n',BBmin(1),BBmin(2),BBmin(3));
    fprintf(fid,'\t%f %f %f\n',BBmax(1),BBmin(2),BBmax(3));
    fprintf(fid,'\t%f %f %f\n',BBmax(1),BBmax(2),BBmax(3));
    fprintf(fid,'\t%f %f %f\n',BBmin(1),BBmax(2),BBmax(3));
    fprintf(fid,'\t%f %f %f\n',BBmin(1),BBmin(2),BBmax(3));

    fprintf(fid,['NBranchingPoints 0\n' ...
                    'NVerticesOnCurves 0\n' ...
                    'BoundaryCurves 0\n' ...
                    'Patches 1\n' ...
                    '{\n' ...
                    'InnerRegion\n' ...
                    'OuterRegion Exterior\n' ...
                    'BoundaryID 0\n' ...
                    'BranchingPoints 0\n\n' ...
                    'Triangles 12\n' ...
                    '  1 6 5\n' ...
                    '  1 2 6\n' ...
                    '  2 7 6\n' ...
                    '  2 3 7\n' ...
                    '  3 8 7\n' ...
                    '  3 4 8\n' ...
                    '  4 1 8\n' ...
                    '  1 5 8\n' ...
                    '  4 3 2\n' ...
                    '  4 2 1\n' ...
                    '  5 6 7\n' ...
                    '  5 7 8\n' ...
                    '}\n']); 
    fclose(fid);
end