%% Return MotifSpectrum of possible triplet motifs
% there are 64 unique motifs, that can be summarized into 16 classes
% Output:
% - MotifSpectrum.Class [64 x 1] (values between 1 and 16)
% - MotifSpectrum.Connectivity [64 x 9] connectivity of triplet motif
%       0-0,0-1,0-2,1-0,1-1,1-2,2-0,2-1,2-2,
% - MotifSpectrum.Label [1 x 16] String with Motif Name 
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
function MotifSpectrum = getMotifSpectrum(filename)    
    if nargin<1
        filename = ['D:\udvary-et-al-2022\analysis\' ...
                        'data\motifs.csv'];
    end
    
    formatSpec = '%d%d%d%d%d%d%d%d%d%d%d%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID,formatSpec,'Delimiter',',', ...
                'HeaderLines',1,'ReturnOnError',false,'EndOfLine','\r\n');
    fclose(fileID);

    connectivity = cell2mat(dataArray(3:11));
    
    %MotifSpectrum.Index = dataArray{1}+1;
    MotifSpectrum.Class = dataArray{2}+1;
    MotifSpectrum.Connectivity = nan(size(connectivity,1),3,3); 
    MotifSpectrum.Label = {'recurrentLoop','directedLoop', ...
                'recurrentIncompleteLoop','directedRecurrentLoop', ...
                'recurrentFFConvergent','recurrentFFDivergent', ...
                'FF','recurrentIncomplete','FFIncomplete', ...
                'recurrentDivergent','recurrentConvergent', ...
                'FFConvergent','FFDivergent', ...
                'recurrentSparse','FFSparse','empty'}; 
    
    for m = 1:size(connectivity,1)
        
        c = 1; 
        
        for i = 1:3
            for j = 1:3
                 MotifSpectrum.Connectivity(m,i,j) = connectivity(m,c); 
                 c = c+1; 
            end
        end
    end
    
end