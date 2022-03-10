% Saves all triplet motif distributions as a table in one mat file
% Outputs: TripletMotifs_All.mat
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
inputPath = [matlabPath 'preprocessing\data\motifs\']; 
outputPath = [matlabPath 'data\']; 

%%
tbl = table();

for i = 1:5

    switch i
        case 1
            strFolder = 'all-column-combinations';
            strPostfix = strFolder;
        case 2
            strFolder = 'selected-column-combinations';
            strPostfix = strFolder;
        case 3
            strFolder = 'combinations_celltype'; 
            strPostfix = 'celltype-combinations';
        case 4
            strFolder = 'intersomatic-distance-combinations';
            strPostfix = strFolder;
        case 5
            strFolder = 'celltype-layer-combinations';
            strPostfix = strFolder;
    end

    tbl_tmp = readtable([inputPath strFolder '\probabilities_16_' strPostfix '.csv']);
    tbl_tmp.pMotif = [tbl_tmp.motif_1_model tbl_tmp.motif_2_model ...
                  tbl_tmp.motif_3_model tbl_tmp.motif_4_model ...
                  tbl_tmp.motif_5_model tbl_tmp.motif_6_model ...
                  tbl_tmp.motif_7_model tbl_tmp.motif_8_model ...
                  tbl_tmp.motif_9_model tbl_tmp.motif_10_model ...
                  tbl_tmp.motif_11_model tbl_tmp.motif_12_model ...
                  tbl_tmp.motif_13_model tbl_tmp.motif_14_model ...
                  tbl_tmp.motif_15_model tbl_tmp.motif_16_model];
    tbl_tmp.pMotifRnd = [tbl_tmp.motif_1_random tbl_tmp.motif_2_random ...
                  tbl_tmp.motif_3_random tbl_tmp.motif_4_random ...
                  tbl_tmp.motif_5_random tbl_tmp.motif_6_random ...
                  tbl_tmp.motif_7_random tbl_tmp.motif_8_random ...
                  tbl_tmp.motif_9_random tbl_tmp.motif_10_random ...
                  tbl_tmp.motif_11_random tbl_tmp.motif_12_random ...
                  tbl_tmp.motif_13_random tbl_tmp.motif_14_random ...
                  tbl_tmp.motif_15_random tbl_tmp.motif_16_random];
    tbl_tmp(:,4:35) = [];

    tblFeatures = readtable([inputPath strFolder '\features_' strPostfix '.csv']);

    % Merge with tblFeatures abg_all and sd_all
    idx = strcmp(tbl_tmp.A,tblFeatures.A) & strcmp(tbl_tmp.B,tblFeatures.B) & ...
                strcmp(tbl_tmp.C,tblFeatures.C);
    if sum(idx)~=length(tblFeatures.A)
        error('Triplet Combinations do not match!');
    end

    tbl_tmp.avg_all = tblFeatures.avg_all;
    tbl_tmp.sd_all = tblFeatures.sd_all; 

    tbl = vertcat(tbl,tbl_tmp); 
end

save([outputPath 'TripletMotifs_All.mat'],'tbl'); 