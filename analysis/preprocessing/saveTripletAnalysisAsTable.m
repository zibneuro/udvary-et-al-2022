% Saves all triplet motif distributions within C2 as a table 
% and adds inDegree-correlation values of respective cells
% Outputs: TripletMotifs_C2.mat
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
dataPath = [matlabPath 'data\']; 
filepath = [matlabPath 'preprocessing\data\motifs\combinations_celltype\'];
tblFeatures = readtable([filepath 'features_celltype-combinations.csv']);  
tblMotifs = readtable([filepath 'probabilities_16_celltype-combinations.csv']);  

%%
tbl = table();
tbl.A = tblMotifs.A;
tbl.B = tblMotifs.B;
tbl.C = tblMotifs.C; 
tbl.pMotif = [tblMotifs.motif_1_model tblMotifs.motif_2_model ...
              tblMotifs.motif_3_model tblMotifs.motif_4_model ...
              tblMotifs.motif_5_model tblMotifs.motif_6_model ...
              tblMotifs.motif_7_model tblMotifs.motif_8_model ...
              tblMotifs.motif_9_model tblMotifs.motif_10_model ...
              tblMotifs.motif_11_model tblMotifs.motif_12_model ...
              tblMotifs.motif_13_model tblMotifs.motif_14_model ...
              tblMotifs.motif_15_model tblMotifs.motif_16_model];
tbl.pMotifRnd = [tblMotifs.motif_1_random tblMotifs.motif_2_random ...
              tblMotifs.motif_3_random tblMotifs.motif_4_random ...
              tblMotifs.motif_5_random tblMotifs.motif_6_random ...
              tblMotifs.motif_7_random tblMotifs.motif_8_random ...
              tblMotifs.motif_9_random tblMotifs.motif_10_random ...
              tblMotifs.motif_11_random tblMotifs.motif_12_random ...
              tblMotifs.motif_13_random tblMotifs.motif_14_random ...
              tblMotifs.motif_15_random tblMotifs.motif_16_random];
          
% Merge with tblFeatures abg_all and sd_all
idx = strcmp(tbl.A,tblFeatures.A) & strcmp(tbl.B,tblFeatures.B) & ...
            strcmp(tbl.C,tblFeatures.C);
if sum(idx)~=length(tblFeatures.A)
    error('Triplet Combinations do not match!');
end

% Clean up cell names
tbl.A = strrep(tbl.A,'C2-','');
tbl.B = strrep(tbl.B,'C2-','');
tbl.C = strrep(tbl.C,'C2-','');

tbl.avg_all = tblFeatures.avg_all;
tbl.sd_all = tblFeatures.sd_all; 

% Add correlations
R = load([dataPath 'CorrData.mat']);
tbl.R_avg = nan(size(tbl.A));
tbl.R = nan(size(tbl.A,1),9); 
for i = 1:size(tbl,1)
    
    % Find corresponding correlation values
    uniCT = unique({tbl.A{i},tbl.B{i},tbl.C{i}}); 
    idx = false(size(R.PreType1)); 
    
    for idx1 = 1:length(uniCT)
        for idx2 = 1:length(uniCT)
            if idx1<=idx2
                continue;
            end
            for idx3 = 1:length(uniCT)
                idx = idx | (strcmp(uniCT{idx1},R.PreType1) & ...
                        strcmp(uniCT{idx2},R.PreType2) & ...
                        strcmp(uniCT{idx3},R.PostType)); 
            end
        end
    end
   
    if sum(idx)==0
       fprintf('%d: %s %s %s\n',sum(idx),tbl.A{i},tbl.B{i},tbl.C{i});  
       continue; 
    end
    
    tbl.R_avg(i) = mean(R.R(idx)); 
    tbl.R(i,1:sum(idx)) = R.R(idx);
end

idxVPM = strcmp(tbl.A,'VPM') & strcmp(tbl.B,'VPM') & strcmp(tbl.C,'VPM');
tbl(idxVPM,:) = []; 

save([dataPath 'TripletMotifs_C2.mat'],'tbl'); 