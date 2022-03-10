function [freqContacts,ratioRange,meanContacts] = ...
                    computeNumSynapses(synapseList,numContacts,ratio,w)
%% Compute the Number of Synpases Per Connections given single cell innervation
% Input:
% - synapseList: vector containing number of innervation between individual
%   cells
% - numContacts (optional): range of contacts
%       default: 1:20
% - ratio (optional): ratio of given range of synaptic contacts 
%       default: 0.99
% - w (optional): weights or frequency for each value in synapseList
%       default: ones(size(synapseList));
% Output:
% - freqContacts: normalized frequency of numContacts
% - ratioRange: range of synapses within defined ratio of 
%   Poission Distribution
% - meanContacts: average number of contacts
%
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)

    if nargin<2
       numContacts = 1:20;
    end
    if nargin<3
        ratio = 0.99;
    end
    if isempty(ratio)
       ratio = 0.99; 
    end
        
    if size(synapseList,1)>1
        synapseList = synapseList';
    end
    
    if nargin<4
       w = ones(size(synapseList));          
    end
    
    if length(w)~=length(synapseList)
       error('w is not same size as synapseList'); 
    end
    
    freqContacts = zeros(size(numContacts)); 
    
    for i = 1:length(synapseList)
        freqContacts = freqContacts + ...
                        w(i).*poisspdf(numContacts,synapseList(i));
    end
    
    % Poisson distribution to get number of contacts per connection
    freqContacts = freqContacts./sum(freqContacts); 
    
    if nargout>1
        [~,ratioRange] = max(freqContacts); 
        while (sum(freqContacts(ratioRange))<ratio*sum(freqContacts))
           ratioRange = [min(ratioRange)-1:1:max(ratioRange)+1];
           ratioRange(ratioRange<1) = [];
           ratioRange(ratioRange>length(freqContacts)) = []; 
           if length(ratioRange) == length(freqContacts)
               break;
           end
        end
    end
    
    if nargout>2
       meanContacts = freqContacts*numContacts'; 
    end
    
%     fprintf('Number of Contacts: M = %.2f [%d %d]\n', ...
%             freqContacts*numContacts',numContacts(ratioRange(1)), ...
%             numContacts(ratioRange(end)));
    
end