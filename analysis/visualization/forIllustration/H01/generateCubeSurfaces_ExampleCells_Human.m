%% Generate 3D Surfaces of example cells of H01
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc
warning on;
matlabPath = 'D:\udvary-et-al-2022\analysis\';
addpath(genpath([matlabPath 'functions\']));
dataPath = [matlabPath 'data\H01\dsc\']; 
cells = [5699471417 6513769803 33330054139];
cubesz = 25000;

%% Write cubes
for i = cells
    [m, ~, ~, BB] = readAmira([dataPath 'all-' ...
            num2str(i) '_binary.am']);
    [~,~,~] = mkdir([dataPath 'Surfaces\']);
    [~,~,~] = mkdir([dataPath 'Surfaces\' num2str(i) '\']);
    writeSurfacesBasedOnDensity(m>0,[dataPath 'Surfaces\' num2str(i) '\'], ...
                BB,cubesz);
    d = dir([dataPath 'Surfaces\' num2str(i) '\*.surf']);
    filenames = strcat([d(1).folder '\'],{d(:).name});
    pt = mergeSurface(filenames,[dataPath 'Surfaces\' num2str(i) '_merged.surf'],cubesz);
    writeSpatialGraphContour(pt,[dataPath 'Surfaces\' num2str(i) '_Box.am']);                      
end