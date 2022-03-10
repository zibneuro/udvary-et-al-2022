%% Predicted connection probability between example cells from H01
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc
warning on;
matlabPath = 'D:\udvary-et-al-2022\analysis\';
addpath(genpath([matlabPath 'functions']));
dataPath = [matlabPath 'data\H01\dsc\']; 
cells = [5699471417 6513769803 33330054139];
strName = 'kji';

%%
for preID = cells
    for postID = cells
        
        if preID==postID
            continue;
        end

        str_ct = [num2str(preID) '-' num2str(postID)];
        filename = [dataPath str_ct '.am']; 
        
        if ~isfile(filename)
           fprintf('%s->%s: 0%%\n',strName(preID==cells), ...
                                        strName(postID==cells));
           continue; 
        end
        
        [I, ~, ~, BB] = readAmira(filename);
        p = 1-exp(-sum(I(:)));
        
        %fprintf('%d->%d: %.1f%%\n',preID,postID,p*100);
        fprintf('%s->%s: %.1f%%\n',strName(preID==cells), ...
                                        strName(postID==cells),p*100);
    end
end