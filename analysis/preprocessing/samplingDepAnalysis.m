%% Calculates connectivity statistics for each cell type combination
% Outputs: .mat files in dataParam folder
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
dataPath = [matlabPath 'data\']; 

%% Initalize variables
PreCellTypes = generateCellTypeList('v2');
PostCellTypes = PreCellTypes(1:end-1);

% For fitting
binsz = 0.01; 
XSample = [0:binsz:1+binsz];
XCount = [0:101];
threshold = 100000; 

%%
c = 1;
for preCTID = 1:length(PreCellTypes)
    for postCTID = 1:length(PostCellTypes)

        preType = PreCellTypes{preCTID};           
        postType = PostCellTypes{postCTID};
        
        load([dataPath 'CellMatrix_C2\' preType '_' postType '.mat']);
            
        % Extract Innervation and probability values
        iValues = I.I(:);
        pValues = 1-exp(-iValues);         
        
        % Compute Parameters
        param_I = computeParameters(iValues);
        param_p = computeParameters(pValues);
                            
        save([dataPath 'dataParam/' ...
            preType '_' postType '.mat'],'param_I','param_p');
        
        fprintf('%d/%d\n',c,length(PreCellTypes)*length(PostCellTypes));
        c = c+1;
    end
end

fclose('all'); 