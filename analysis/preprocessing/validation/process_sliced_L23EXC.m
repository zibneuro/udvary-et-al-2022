%% Read in raw data and save 20 matrices (one per slice)
% Generates .mat files for each slice with L23 EXC neurons
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

%% L23EXC-L23EXC (according to Averman et al)
matlabPath = 'D:\udvary-et-al-2022\analysis\';
outputPath = [matlabPath 'preprocessing\data\L23EXC_L23EXC_slice\']; 
id = '3410f9b8-f11a-49c2-b033-f63eb98b2dc9';
inputPath = [matlabPath 'preprocessing\data\dsc_L23\'];
        
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
            
    save([outputPath 'slice-' num2str(sliceID) '.mat'],'I');
    fprintf('Slice %d done\n',sliceID);
end