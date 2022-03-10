%% Compute how many axons/dendrites come from cells whose somata is 
% inside cube, 
% outside cube (but within vS1), 
% outside cube (but from VPM),
% outside cube (outside of vS1 and not in VPM),
% per cube (100x100x50). Values will be different for different volumes
% -> higher percentage when looking at entire volume (see
% outsideInnervation.m)
% Save results in csv files
%
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
addpath([matlabPath 'functions\']);
tablesPath = [matlabPath 'output\Tables\']; 
post = load([matlabPath 'data\dataSubvolume\post'  ...
            '_CellsInsideOut_100_100_50.mat'],'cellInsideOut');
pre = load([matlabPath 'data\dataSubvolume\pre'  ...
            '_CellsInsideOut_100_100_50.mat'],'cellInsideOut');
numCubes = size(post.cellInsideOut,1);

%% Save raw dat as table [pre and post]
tbl = table([1:numCubes]',post.cellInsideOut(:,1),post.cellInsideOut(:,2), ...
            post.cellInsideOut(:,4), ...
            pre.cellInsideOut(:,1),pre.cellInsideOut(:,2), ...
            pre.cellInsideOut(:,3),pre.cellInsideOut(:,4), ...
            'VariableNames',{'Cube','dendInsideCube','dendOutsideCube_vS1', ...
            'dendOutsideCube_rest', ...
            'axonInsideCube','axonOutsideCube_vS1', ...
            'axonOutsideCube_VPM','axonOutsideCube_rest'});
writetable(tbl,[tablesPath 'CubeStats_InsideOut_100x100x50_Raw.csv']);

%% Save summary stats
% Delete outside vS1 innervation (not VPM)
post.cellInsideOut(:,4) = [];
pre.cellInsideOut(:,4) = [];
postVal = bsxfun(@rdivide,post.cellInsideOut,sum(post.cellInsideOut,2));
preVal = bsxfun(@rdivide,pre.cellInsideOut,sum(pre.cellInsideOut,2));
total = pre.cellInsideOut + post.cellInsideOut;
totalVal = bsxfun(@rdivide,total,sum(total,2));

if sum(postVal(:,3))>0
    error('something is weird, there should be no VPM wrt to dendrites');
end

% -> mean median prctl min max
optStringParam = {'dendInsideCube[%]','dendOutsideCube_vS1[%]', ...
                'axonInsideCube[%]','axonOutsideCube_vS1[%]', ...
                'axonOutsideCube_VPM[%]', ...
                'totalInsideCube[%]','totalOutsideCube_vS1[%]', ...
                'totalOutsideCube_VPM[%]'};

fid = fopen([tablesPath 'CubeStats_InsideOut_100x100x50_Summary.csv'],'w');
if fid==-1
    error('Cannot write file');
end
fprintf(fid,'Param,Mean,SD,Median,25th,75th,min,max,n,\n');

for optParam = 1:length(optStringParam)

    switch optParam
        case 1
           x = postVal(:,1)*100; % (dend) inside
        case 2
           x = postVal(:,2)*100; % (dend) outside vS1
        case 3
           x = preVal(:,1)*100; % (axon) inside 
        case 4
           x = preVal(:,2)*100; % (axon)outside vS1
        case 5
           x = preVal(:,3)*100; % (axon)outside VPM
        case 6
           x = totalVal(:,1)*100; % (all) inside    
        case 7
           x = totalVal(:,2)*100; % (all) outside vS1
        case 8
           x = totalVal(:,3)*100; % (all) outside VPM
    end

    fprintf(fid,'%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,\n', ...
                optStringParam{optParam},mean(x),std(x), ...
                median(x),prctile(x,25),prctile(x,75), ...
                min(x),max(x),numel(x));
end
fclose('all');

optParam = 8;
x = totalVal(:,1)*100; % (all) inside
fprintf('Param,Mean,SD,Median,25th,75th,min,max,n,\n');
fprintf('%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,\n', ...
            optStringParam{optParam},mean(x),std(x), ...
            median(x),prctile(x,25),prctile(x,75), ...
            min(x),max(x),numel(x));