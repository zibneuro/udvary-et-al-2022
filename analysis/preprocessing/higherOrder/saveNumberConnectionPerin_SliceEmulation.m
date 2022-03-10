%% Saves higher order connectivity patterns for comparison with Perin et al., 2011
% Number of connections in clusters (Fig 2A and 2B of Perin et al., 2011)
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all;
close all;
clc;

matlabPath = 'D:\udvary-et-al-2022\analysis\';
inputPath = [matlabPath 'preprocessing\data\L5PT_L5PT_slice\']; 
outputPath = [matlabPath 'data\higherOrder\dataPerin\']; 
sliceIDs = [0:1:19];
numTrials = 1000;          % number of neurons extracted
samplingThreshold = 1000;   % number of random samples for each edge 
                            % if full spectrum is not computed 

%%
rng(100842);                                
for k = 1:length(sliceIDs)

    load([inputPath 'slice-' num2str(sliceIDs(k)) '.mat']);
    pMatrix = 1-exp(-I.I);
    n = size(pMatrix,1); 
    p_avg = mean(pMatrix(:)); 
 
    for numNodes = 3:8 

        numDirEdges = numNodes * numNodes - numNodes; 
        strNumNodes = num2str(numNodes); 

        p = nan(numTrials,numNodes,numNodes);

        % Matrix Mapper
        id = 1:(numNodes*numNodes); % all ids
        idDiagonal = 1:numNodes+1:numNodes*numNodes; % id of diagonal values
        idValid = uint8(setdiff(id,idDiagonal)); 

        for j = 1:numTrials
            idx = randperm(n,numNodes); 
            ptmp = pMatrix(idx,:);
            ptmp = ptmp(:,idx);
            ptmp(idDiagonal) = nan; 
            p(j,:,:) = ptmp; 
        end

        p_avg_sample = nanmean(p(:)); 

        pMotifUniform = nan(1,numDirEdges+1); 
        pMotif = zeros(numTrials,numDirEdges+1); 
        numCombinationsList = nan(1,numDirEdges/2+1);

        fprintf('------------------\n');
        fprintf('>> %d Cells (%d directed Edges) (slice %d)\n', ...
                numNodes, numDirEdges, sliceIDs(k));

        for currNumDirEdges = 0:numDirEdges/2

            numCombinations = nchoosek(numDirEdges,currNumDirEdges); 
            numCombinationsList(currNumDirEdges+1) = numCombinations;
            currNumDirEdgesCompl = numDirEdges-currNumDirEdges; 

            fprintf('currNumEdges %d & %d: #Combinations = %d', ...
                        currNumDirEdges,currNumDirEdgesCompl,numCombinations);

            if numCombinations>samplingThreshold

                fprintf(' >> Sampling (%d)\n',samplingThreshold);

                pMotif = computeSampling(numDirEdges,currNumDirEdges, ...
                                        samplingThreshold,idValid,numNodes, ...
                                        numTrials,p,pMotif);
                normConstant = 1; % not entire spectrum, only average sample

                % If multiplying setting normConstant to actual numCombinations
                % here, then pMotifUniform will sum up to 1
                % however, pMotif model will not, so one has to implement that
                % pMotif is also somehow scaled by normConstant (since there 
                % are only samples, it has to be done in computeSampling)
            else
                fprintf('\n');
                pMotif = computeCombinations(numDirEdges,currNumDirEdges, ...
                                    numCombinations,idValid,numNodes, ...
                                    numTrials,p,pMotif);
                normConstant = numCombinations; % entire spectrum
            end

            pMotifUniform(currNumDirEdges+1) = p_avg^currNumDirEdges * ... 
                        (1-p_avg)^(numDirEdges-currNumDirEdges) * normConstant; 

            if currNumDirEdges ~= numDirEdges/2
                pMotifUniform(currNumDirEdgesCompl+1) = p_avg^currNumDirEdgesCompl * ... 
                        (1-p_avg)^(numDirEdges-currNumDirEdgesCompl) * normConstant;  
            end
        end

        pMotif_avg = mean(pMotif); 
        pDev = (pMotif_avg-pMotifUniform)./pMotifUniform; % Perin Deviation

        % Check whether sum is equal to 1, if not, something is wrong
        % However, can be due to sampling
        if abs(sum(pMotif_avg)-1)>1e-10
            warning('Sum of probabilities is not equal to 1 but %.2e', ...
                        sum(pMotif_avg));
        end

        if abs(sum(pMotifUniform)-1)>1e-10
            warning('Sum of expected probabilities is not equal to 1 but %.2e', ...
                        sum(pMotifUniform));
        end

        save([outputPath strNumNodes '_' num2str(numTrials) '_Sampling_' ...
                num2str(samplingThreshold) '_Slice-' ...
                num2str(sliceIDs(k)) '.mat'], 'pMotif_avg', ...
                'pMotifUniform','pDev','p_avg_sample','numDirEdges', ...
                'numNodes','numTrials','samplingThreshold', ...
                'numCombinationsList'); 
    end
end

%% Functions
function pMotif = computeCombinations(numDirEdges,currNumDirEdges, ...
                                    numCombinations,idValid,numNodes, ...
                                    numTrials,p,pMotif)

    ii = nchoosek(uint8(1:numDirEdges),uint8(currNumDirEdges)); 
    idMatrix = idValid(ii);
    currNumDirEdgesCompl = numDirEdges-currNumDirEdges; 
    
    for i = 1:numCombinations

        m = false(numNodes,numNodes); 

        if ~isempty(idMatrix)
            if size(idMatrix,1)==1
                if numCombinations==1 
                    % here all indices for one combination 
                    % (fully connected)
                    m(idMatrix) = 1; 
                else        
                    % here each indices for one combination 
                    % (fully unconnected)
                    m(idMatrix(i)) = 1; 
                end
            else
                m(idMatrix(i,:)) = 1; 
            end
        end

        if currNumDirEdges ~= numDirEdges/2
            mCompl = ~m; % Complement
        end

        % probability of combination
        for j = 1:numTrials
            ptmp = squeeze(p(j,:,:)); 
            pMotif(j,currNumDirEdges+1) = ...
                            pMotif(j,currNumDirEdges+1) ...
                            + prod(ptmp(m & ~isnan(ptmp))) ...
                            * prod(1-ptmp(~m & ~isnan(ptmp)));

            if currNumDirEdges ~= numDirEdges/2
                pMotif(j,currNumDirEdgesCompl+1) = ...
                                pMotif(j,currNumDirEdgesCompl+1) ...
                                + prod(ptmp(mCompl & ~isnan(ptmp))) ...
                                * prod(1-ptmp(~mCompl & ~isnan(ptmp)));
            end
        end
    end
end


function pMotif = computeSampling(numDirEdges,currNumDirEdges, ...
                                    numSamples,idValid,numNodes, ...
                                    numTrials,p,pMotif)
        
    currNumDirEdgesCompl = numDirEdges-currNumDirEdges; 
                                
    for i = 1:numSamples

        m = false(numNodes,numNodes); 

        ii = randperm(length(idValid),currNumDirEdges); 
        idMatrix = idValid(ii);

        if ~isempty(idMatrix)
            m(idMatrix) = 1; 
        end

        if currNumDirEdges ~= numDirEdges/2
            mCompl = ~m; % Complement
        end

        % probability of combination
        for j = 1:numTrials
            ptmp = squeeze(p(j,:,:)); 
            pMotif(j,currNumDirEdges+1) = ...
                            pMotif(j,currNumDirEdges+1) ...
                            + prod(ptmp(m & ~isnan(ptmp))) ...
                            * prod(1-ptmp(~m & ~isnan(ptmp)));

            if currNumDirEdges ~= numDirEdges/2
                pMotif(j,currNumDirEdgesCompl+1) = ...
                                pMotif(j,currNumDirEdgesCompl+1) ...
                                + prod(ptmp(mCompl & ~isnan(ptmp))) ...
                                * prod(1-ptmp(~mCompl & ~isnan(ptmp)));
            end
        end
    end
    
    pMotif(:,currNumDirEdges+1) = ...
                pMotif(:,currNumDirEdges+1)./numSamples;
    
   if currNumDirEdges ~= numDirEdges/2
        pMotif(:,currNumDirEdgesCompl+1) = ...
                    pMotif(:,currNumDirEdgesCompl+1)./numSamples;
   end
end