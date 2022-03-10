function pMotif = computeMotifProbabilitiesBasedOnDoubletMotifs(pDoublet)
%% Compute probability of each motif given a doublet motif probabilities
% NOTE: unidirectional motif for one direction! (so divide by 2)
% sum should not be 1.
% Input:
% - pDoublet: 1x3
%   pDoublet(1): unidirectinal (for one direction!!!)
%   pDoublet(2): bidirectional
%   pDoublet(3): unconnected
% - pMotif: 1 x 16 with probability of each motif occuring
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)

    if length(pDoublet)~=3
       error('requires 3 values!'); 
    end

    if (min(pDoublet)<0 || min(pDoublet)>1)
        error('probability has to be between 0 and 1')
    end
    
    pMotif = nan(1,16);
    pMotif(1) = pDoublet(2)^3;
    pMotif(2) = pDoublet(2)^2*pDoublet(1);
    pMotif(3) = pDoublet(2)*pDoublet(1)^2;
    pMotif(4) = pDoublet(2)*pDoublet(1)^2;
    pMotif(5) = pDoublet(2)*pDoublet(1)^2;
    pMotif(6) = pDoublet(1)^3;
    pMotif(7) = pDoublet(1)^3;
    pMotif(8) = pDoublet(2)^2*pDoublet(3);
    pMotif(9) = pDoublet(2)*pDoublet(1)*pDoublet(3);
    pMotif(10) = pDoublet(2)*pDoublet(1)*pDoublet(3);
    pMotif(11) = pDoublet(1)^2*pDoublet(3);
    pMotif(12) = pDoublet(1)^2*pDoublet(3);
    pMotif(13) = pDoublet(1)^2*pDoublet(3);
    pMotif(14) = pDoublet(2)*pDoublet(3)^2;
    pMotif(15) = pDoublet(1)*pDoublet(3)^2; 
    pMotif(16) = pDoublet(3)^3; 
    
    % number of permutation of each motif (number of possibilities)
    numPerm = [1,6,3,3,6,6,2,3,6,6,3,3,6,3,6,1]; 
    pMotif = pMotif.*numPerm;

    if abs(sum(pMotif)-1)>50*eps
       warning(['Sum not equal to 1 (%.2e)! Mabye pDoublet was not ' ...
           'normalized correctly! pDoublet(1) for one direction! ' ...
           'Try dividing by two!'],sum(pMotif)); 
    end
end