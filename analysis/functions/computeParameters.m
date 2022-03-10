function [param] = computeParameters(x,modeNumDigits)
%% [param] = computeParameters(x)
% Compute parameters for x (innervation or probability values)
% Input:
% - x: vector [1 x N] containing either innervation or probability values
% - modeNumDigits (optional): Digits after comma used to compute mode
%                   default: 4
% Output:
% - param: structure 
%       SummaryStatistics
%       - param.m: mean
%       - param.SD: standard deviation
%       - param.med: median
%       - param.mod: mode (rounded to 4 digitis after comma)
%       - param.skew: (mean-mode)/SD
%       - param.min: minimum value
%       - param.max: maximum value
%       - param.prctl: percentile distribution [1 5 10 25 75 90 95 99]%
%       - param.n: Number of samples
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)

    if ~isvector(x)
        x = x(:);
    end
    
    if size(x,2)>1
        x = x'; % only works for column vectors
    end
    
    if nargin<2
       modeNumDigits = 4; 
    end
    
    % Mean+-SD,Skewness,SAE(best),SAE(Gauss)
    param.m = mean(x);
    param.SD = std(x);
    param.med = median(x);
    param.mod = mode(round(x,modeNumDigits)); 
    param.skew = (param.m - param.mod)./param.SD;
    param.min = min(x);
    param.max = max(x);
    param.prctl = [prctile(x,1) prctile(x,5) prctile(x,10) prctile(x,25) ...
                    prctile(x,75) prctile(x,90) prctile(x,95) prctile(x,99)];
    param.n = numel(x);    
end