function [motifProb16] = collapse64MotifTo16Motif(motifProb64)
%% [motifProb16] = collapse64MotifTo16Motif(motifProb64)
% Transforms 64 motif probability in 16 motif probability by summing 
% up probabilities
% Input:
%   - motifProb64 [64x1 or 1x64] probabilities of each motif
%       or  [64xN] probabilities of each motif for N samples
% Output:
%   - motifProb16 [16x1] probability of each summarized motif
%       or motifProb16 [16xN] probability of each summarized motif
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)

    numMotifs = 16;
    idx = [1 2 2 3 3 3 3 3 3 4 4 4 4 4 4 5 5 5 6 6 6 7 7 7 7 7 7 8 8 8 ...
        9 9 9 9 9 9 10 10 10 10 10 10 11 11 11 11 11 11 12 12 12 ...
        13 13 13 14 14 14 15 15 15 15 15 15 16]; 

    if numel(motifProb64)==64
    
        motifProb16 = nan(numMotifs,1); 

        for i = 1:numMotifs
            motifProb16(i) = sum(motifProb64(idx==i)); 
        end

    elseif size(motifProb64,1)==64 && size(motifProb64,2)>1 
        
        numSamples = size(motifProb64,2);
        motifProb16 = nan(numMotifs,numSamples); 

        for i = 1:numMotifs
            motifProb16(i,:) = sum(motifProb64(idx==i,:),1); 
        end
        
    else
       error('motifProb64 has be [64x1, 1x64, or 64xN]!'); 
    end
    
end