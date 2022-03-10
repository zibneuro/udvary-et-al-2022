%% Read in raw data and save 20 matrices (one per slice)
% Generates .mat files for each slice with L5PT neurons
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

%% L5PT-L5PT (according to Perin et al)
matlabPath = 'D:\udvary-et-al-2022\analysis\';
outputPath = [matlabPath 'preprocessing\data\L5PT_L5PT_slice\']; 
id = '8cb7701b-ed14-4743-9ba3-d4a84d3b585d';
inputPath = [matlabPath 'preprocessing\data\dsc_L5PT\'];

        
%%
for sliceID = 0:19 

    d = dir([inputPath id '_slice-' num2str(sliceID) '\*_DSC.csv']);

    % Extract values
    tbl = table(); 
    for i = 1:length(d)
        tbltmp = readtable([d(i).folder '\' d(i).name],'Format','%d%f');
        pre_id = str2double(regexp(d(i).name,'\d+','match'));
        tbltmp.pre_id = pre_id.*int32(ones(size(tbltmp.post_id))); 
        tbl = [tbl; tbltmp]; 
    end

    % Create Matrix
    I.PreCellID = unique(tbl.pre_id)';
    I.PostCellID = unique(tbl.post_id);
    I.I = zeros(length(I.PreCellID),length(I.PostCellID)); 

    [~,idxPre] = ismember(tbl.pre_id,I.PreCellID);
    [~,idxPost] = ismember(tbl.post_id,I.PostCellID);
    idx = sub2ind(size(I.I),idxPre,idxPost);
    I.I(idx) = tbl.dsc; 
            
    save([outputPath '\slice-' num2str(sliceID) '.mat'],'I');
    fprintf('Slice %d done\n',sliceID);
end