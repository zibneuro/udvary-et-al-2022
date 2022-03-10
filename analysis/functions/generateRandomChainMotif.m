function edgeMatrix = generateRandomChainMotif(numNodes,numSamples)
% edgeMatrix = generateRandomChainMotif(numNodes,numSamples)
% Generates binary edge matrices with chain motifs
% Input:
% - numNodes: number of nodes (size of square matrix)
% - numSamples: number of motifs/samples
% Output:
% - edgeMatrix: binary matrix [numSamples x numNodes x numNodes]
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)

    numDirEdges = numNodes-1;
    
    if nargin==1
        numSamples = 1;
    end

    edgeMatrix = false(numSamples,numNodes,numNodes);
    
    for t = 1:numSamples
        idxValidRows = 1:numNodes; 
        idxValidCol = 1:numNodes; 
        currEdgeMatrix = false(numNodes,numNodes); 
        
        for i = 1:numDirEdges

            r = idxValidRows(randperm(numel(idxValidRows),1)); 

            idxtmp = idxValidCol;
            idxtmp(idxtmp==r) = []; 

            c = idxtmp(randperm(numel(idxtmp),1)); 

            currEdgeMatrix(c,r) = true; 

            idxValidRows(idxValidRows==r) = [];
            idxValidCol(idxValidCol==c) = [];
        end

        % Double checking
        if (numDirEdges~=sum(sum(currEdgeMatrix)==1) || ...
                numDirEdges~=sum(sum(currEdgeMatrix,2)==1) || ...
                1~=sum(sum(currEdgeMatrix)==0) || ...
                1~=sum(sum(currEdgeMatrix,2)==0))
           error('Current edge motif is not a chain motif!');  
        end
        
        edgeMatrix(t,:,:) = currEdgeMatrix;
    end
end