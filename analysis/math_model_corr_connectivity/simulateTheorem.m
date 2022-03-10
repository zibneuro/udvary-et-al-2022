function [result] = simulateTheorem(numSamples,gamma,lambda,numEdges, ...
                                                p_mu_analytical,edgeList)
%% [result] = simulateTheorem(numSamples,lambda,gamma)
% Input:
% - numSamples: number of Samples drawn from shared Source S
% - gamma: mean of shared source
% - lambda: variance of shared source
% - numEdges: Number of Edges investigated (Motif size / Number of nodes)
% - p_mu_analytical (optional): To test imprecision
% - edgeList (optional): Number of edges to iterate over
%       default: 1:numEdges.
% Output:
% - result.p_mu: mean(p)
% - result.p_var: var(p); 
% - result.pMotif: Motif probability [1 x numEdges];
% - result.pMotif0: Motif probability if numEdges==0 
% - result.pRandom: Motif probability only mean [1 x numEdges];
% - result.pRandom0: Motif probability only mean and if numEdges==0 
% - result.sum_imprecision: Imprecision with respect to sum over
%               probabilities (should be 1)
% - result.mean_imprecision: Imprecision with respect to mean(p) and
%               analytical mean computed as integral
%
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
% Date: Feb, 9th, 2022

    if nargin < 5
        p_mu_analytical = [];
    end
    if nargin<6
       edgeList = 1:numEdges; 
    end
    if max(edgeList)>numEdges
       error('Maximum in edgeList (%d) is larger than numEdges (%d)', ...
           max(edgeList),numEdges); 
    end
    if min(edgeList)<1
       error('Minimum in edgeList (%d) has to be larger than 0', ...
           min(edgeList)); 
    end

    % Get norm cdf (note as input standard deviations!)
    L = @(x,mu,sd) 1 - normcdf(x,mu,sd);
    
    % Variance of private source T
    n = 1 - lambda; 

    % Shared source S
    S = normrnd(0,1,numSamples,1);         

    % ----------------------------------
    % Derived directly from S
    % Compute connection probability given S
    
    if (lambda==0)
       p = L(0,gamma,sqrt(n));  % 1 x 1
       % to avoid issues with var
       % in this case var(p) = 0
       % and pMotif == pRandom
    else
       p = L(0,gamma+sqrt(lambda).*S,sqrt(n)); % 1 x numSamples
    end

    % Mean over connection probabilities
    p_mu = mean(p); 

    % Variance over connection probabilities
    p_var = var(p); 

    % Motif probabilities       
    pMotif = nan(1,numEdges); 
    pRandom = nan(1,numEdges); 
    
    for k = edgeList
        
        % Binomial Coefficient
        bc = nchoosek(numEdges,k); 

        % Motif probability based on S
        pMotif(k) = bc*mean(p.^k .* (1-p).^(numEdges-k));

        % Motif probability based on mean
        pRandom(k) = bc*(p_mu^k*(1-p_mu)^(numEdges-k));
        
    end

    % Unconnected motif
    pRandom0 = (1-p_mu)^(numEdges);
    pMotif0 = mean((1-p).^numEdges);
        
    % Results    
    if isempty(p_mu_analytical) || isnan(p_mu_analytical)
        p_mu_analytical = L(0,gamma,sqrt(lambda+n));
    end

    result.p_mu = p_mu;
    result.p_var = p_var; 
    result.pMotif = pMotif;
    result.pMotif0 = pMotif0; 
    result.pRandom = pRandom; 
    result.pRandom0 = pRandom0; 
    result.mean_imprecision = abs(p_mu-p_mu_analytical); 
    
    % Only check if iterated over all edges
    if sum(abs(edgeList-[1:numEdges]))==0
        % Check imprecisions
        % Check whether sum is zero 
        %   (therefore compute prob of unconnected motif)
        pRandom_sum_check = sum(pRandom) + pRandom0;
        pMotif_sum_check = sum(pMotif) + pMotif0;

        if sum(abs(pRandom_sum_check-1)>100*eps)>0
            warning(['Sum of pRandom is not equal to 1.0! ' ...
                '(diff = %.2e, gamma = %.2f, lambda = %.2f)'], ...
                abs(pRandom_sum_check-1),gamma,lambda);
        end

        if abs(pMotif_sum_check-1)>100*eps
            warning(['Sum of pMotif is not equal to 1.0! ' ...
                '(diff = %.2e, gamma = %.2f, lambda = %.2f)'], ...
                abs(pMotif_sum_check-1),gamma,lambda);
        end
        
        if abs(pMotif_sum_check-1)>min([pMotif pMotif0])/100
            warning(['Sum of pMotif is less precise than min(pMotif) ' ...
                '(%.2e vs. %.2e) (gamma = %.2f, lambda = %.2f)!'], ...
                        abs(pMotif_sum_check-1),min([pMotif pMotif0]), ...
                        gamma,lambda);
        end
        
        if abs(pRandom_sum_check-1)>min([pRandom pRandom0])/100
            warning(['Sum of pRandom is less precise than min(pRandom) ' ...
                '(%.2e vs. %.2e) (gamma = %.2f, lambda = %.2f)!'], ...
                        abs(pRandom_sum_check-1),min([pRandom pRandom0]), ...
                        gamma,lambda);
        end

        result.sum_imprecision = max([abs(pMotif_sum_check-1) ...
                                        abs(pRandom_sum_check-1)]); 
    else
       result.sum_imprecision = nan;  
    end
end