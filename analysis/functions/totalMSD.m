%% Computes Total/Pooled Mean and SD over different Subpopulations
% Input:
%   - m: mean [1xn] (optional; default: [])
%   - s: sd [1xn] (optional; default: []) (unbiased estimator)
%   - n: sample size [1xn]
% Output:
%   - M: pooled Mean (if no m given, M is empty)
%   - SD: pooled SD (if no s and no m given, SD is empty)
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
%
% Keywords: Pooled Standard Deviation/SD/STD, Pooled Mean/Average,
function [M,SD] = totalMSD(m,s,n)

    M = [];
    SD = []; 

    if isempty(n)
        error('Size of Samples n is missing')
    elseif size(n,1)>1
        n = n'; 
    end
    
    if sum(isnan(m))>0
        error('Mean m contains NaN')
    end
    if sum(isnan(s))>0 
        error('Standard Deviation s contains NaN')
    end
    if sum(isnan(n))>0
        error('Sample size n contains NaN')
    end
    
    if numel(unique(n))==1
        warning(['Sample size is identical for each sample!' ...
            ' Do not use this function, take mean instead!']);
    end
    
    if ~isempty(m)
        if length(m) ~= length(n)
            error('Sample size n and mean values m have different lenghts'); 
        end
        if size(m,1)>1
            n = m'; 
        end
        if length(m)==1
            M = m; 
        else
            M = m*n'/sum(n); 
        end
    end
    
    if ~isempty(s) && ~isempty(m)
        if length(s) ~= length(n)
            error('Sample size n and sd values s have different lenghts'); 
        end
        if size(s,1)>1
            n = s'; 
        end
        
        if length(s)==1
            SD = s;
        else
            % Take care of n=0
            s(n<1) = [];
            m(n<1) = [];
            n(n<1) = [];

            m(isnan(s)) = [];
            n(isnan(s)) = [];
            s(isnan(s)) = [];

            % Calculate SD
            s = s .* sqrt((n-1)./n); % convert to biased estimator       
            var_tmp = (n)*(s.^2+m.^2)' / sum(n) - (m*n'/sum(n))^2 ; % calculate variance
            SD = sqrt(var_tmp) * sqrt(sum(n)/(sum(n)-1)); % convert to unbiased estimator 
        end
    end
end