function [infGamma,infLambda,result] = ...
                    getTheoremPrediction(filename,p_avg,p_var)
%% [gammaValue,lambdaValue] = getTheoremPrediction(p_avg,p_var)
% Uses lookup table to get the corresponding/inferred lambda and gamma
% value for any average and variance of connection probabilities
% Input:
% - filename: /path/to/Simulation_dd_mm_yyyy_HH_MM.mat
% - p_avg: average connection probability
% - p_var: variance of connection probability
% Output:
% - infGamma: inferred gamma value
% - infLambda: inferred lambda value
% - result: 
%     result.diffVarValues: Difference between variance in lookup table and
%     actual variance
%     result.diffMeanValues: Difference between mean in lookup table and
%     actual mean
%     result.devMotif: [1 x 7] theorem-based deviation for each edge number
%           in motif of size 3. (0:6)
%
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
% Date: Feb, 9th, 2022

    load(filename,'gamma','var_sample','mu_analytical','lambda','devMotif'); 
    var_sample = mean(var_sample,3);
    
    % mu is identitical across columns
    [~,idxMU] = min(abs(mu_analytical(:,1)-p_avg));

    % find variance given mu
    [~,idxVAR] = min(abs(var_sample(idxMU,:)-p_var));
    
    % Differences between values
    result.diffVarValues = abs(var_sample(idxMU,idxVAR)-p_var);
    result.diffMeanValues = abs(mu_analytical(idxMU,idxVAR)-p_avg);
    result.devMotif = squeeze(devMotif(idxMU,idxVAR,:)); 

    % Get corresponding gamma and lambda values
    infGamma = gamma(idxMU,idxVAR);
    infLambda = lambda(idxMU,idxVAR); 
end