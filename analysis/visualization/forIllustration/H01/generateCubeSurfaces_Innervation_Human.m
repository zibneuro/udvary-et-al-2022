%% Generate 3D Surfaces of the overlap of the example cells of H01
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc
warning on;
matlabPath = 'D:\udvary-et-al-2022\analysis\';
addpath(genpath([matlabPath 'functions\']));
dataPath = [matlabPath 'data\H01\dsc\']; 
surfacePath = [dataPath 'Surfaces\Innervation\'];
[~,~,~] = mkdir(surfacePath);
cells = [5699471417 6513769803 33330054139];
cubesz = 25000;

%%
for preID = cells
    for postID = cells
        
        if preID==postID
            continue;
        end

        str_ct = [num2str(preID) '-' num2str(postID)];
        filename = [dataPath str_ct '.am']; 
        
        if ~isfile(filename)
           continue; 
        end
        
        [I, ~, ~, BB] = readAmira(filename);
        p = 1-exp(-I); 

        % All overlapping cubes!
        [~,~,~] = mkdir([surfacePath str_ct '\']);
        writeSurfacesBasedOnDensity(p>0,[surfacePath str_ct '\'],BB,cubesz);
        d = dir([surfacePath str_ct '\*.surf']);
        filenames = strcat([d(1).folder '\'],{d(:).name});
        pt = mergeSurface(filenames,[surfacePath str_ct '_merged.surf'],cubesz);
        writeSpatialGraphContour(pt,[surfacePath str_ct '_Box.am']);

        % Cubes per value
        numCubes = 5;
        pval = unique(p(:));
        x = linspace(min(pval),max(pval),numCubes);
        fprintf('---- %s ----\n',str_ct);
        fprintf('%d overlapping cubes\n',sum(p(:)>0));
        fprintf('p_max = %.2f%% [%.2e%%]\n',max(pval)*100,max(pval)*100);

        Y = discretize(p,x);
        Y(p==0) = 0;
        cube_list = unique(Y)';

        for y = cube_list(2:end)
            strtmp = sprintf('%s_cube_%d',str_ct,y);
            [~,~,~] = mkdir([surfacePath strtmp '\']);
            writeSurfacesBasedOnDensity(Y==y,[surfacePath strtmp '\'],BB,cubesz);
            d = dir([surfacePath strtmp '\*.surf']);
            filenames = strcat([d(1).folder '\'],{d(:).name});
            pt = mergeSurface(filenames,[surfacePath strtmp '_merged.surf'],cubesz);
            %writeSpatialGraphContour(pt,[surfacePath str '_Box.am']);
        end
    end
end