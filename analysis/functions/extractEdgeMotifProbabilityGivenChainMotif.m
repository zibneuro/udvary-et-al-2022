function [edgeStats] = extractEdgeMotifProbabilityGivenChainMotif( ...
                                edgeMatrix,pSelect,p_avg)  
% [edgeStats] = extractEdgeMotifProbabilityGivenChainMotif( ...
%                                               edgeMatrix,pSelect,p_avg)    
% Compute edgeStats given set of motifs
% Input: 
% - edgeMatrix [numMotifs x numNodes x numNodes]
%       binary matrix with 0s and 1s for connecting edges
% - pSelect [numTrials x numNodes x numNodes]
% - p_avg
% Output:
% - edgeStats
%     edgeStats.numCells 
%     edgeStats.numMotifs 
%     edgeStats.pEdgeUniformMotif
%     edgeStats.pEdgeMotif [numTrials x numMotifs]
%     edgeStats.numOfTrials
%     edgeStats.p_avg
% NOTE: THIS is only relative probability for one chain motif, 
%       not for all combinations that yield a chain motif
%       To get absolute probability of all chain motifs, multiply
%       probability by binomial coefficient
%       bc = nchoosek(numEdges,k); 
%           k = current number of edges, here numNodes-1
%           numEdges = numNodes*numNodes-numNodes (total number of possible
%           edges)
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)

    numMotifs = size(edgeMatrix,1); 
    numTrials = size(pSelect,1);
    numNodes = size(pSelect,2); 

    pMotif = nan(numTrials,numMotifs);    
    connected = numNodes-1;
    unconnected = numNodes*numNodes-numNodes-(numNodes-1);
    % = numNodes * numNodes - 2 * numNodes + 1);
    pMotifUniform = p_avg^connected * (1-p_avg)^unconnected;
    
    for i = 1:numMotifs
        e = squeeze(edgeMatrix(i,:,:));        
        for j = 1:numTrials
            ptmp = squeeze(pSelect(j,:,:)); 
            pMotif(j,i) = prod(ptmp(e & ~isnan(ptmp))) ...
                                * prod(1-ptmp(~e & ~isnan(ptmp)));
        end
    end
    
    edgeStats.numCells = numNodes; 
    edgeStats.numMotifs = numMotifs; 
    edgeStats.pEdgeUniformMotif = pMotifUniform;
    edgeStats.pEdgeMotif = pMotif;
    edgeStats.numOfTrials = numTrials;
    edgeStats.p_avg = p_avg; 
    
end