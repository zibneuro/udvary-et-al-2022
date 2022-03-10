%% generate CellArray that contains all CellTypes
% Input:
% - opt (optional): 
%       '': for all cell types (default) (excl. VPM)
%       inhibitory: all inhibitory cell types (incl LocalSubtypes)
%       excitatory: all excitatory cell types (excl. VPM)
%       ID: all celltypes sorted according to NeuroNet(ZIBAmira)
% Output:
% - CellTypeList: Array containing string cell array of all cell types
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
function [CellTypeList] = generateCellTypeList(opt)

    if nargin==0
        opt = '';
    end
    
    switch opt
        case 'excitatory'
            CellTypeList = {'L2','L34','L4py','L4sp','L4ss','L5st', ...
                'L5tt','L6cc','L6ccinv','L6ct'};
        case 'inhibitory'
            CellTypeList = {'SymLocal1','SymLocal2','SymLocal3', ...
                'SymLocal4','SymLocal5','SymLocal6','L1','L23Trans', ...
                'L45Sym','L45Peak','L56Trans'};
        case 'ID'
            CellTypeList = {'L2','L34','L4py','L4sp','L4ss', ...
                'L5st','L5tt','L6cc','L6ccinv','L6ct', ...
                'VPM','L2axon','L34axon','L4pyaxon','L4spaxon','L4ssaxon', ...
                'L5staxon','L5ttaxon','L6ccaxon','L6ccinvaxon','L6ctaxon', ...
                'SymLocal','SymLocal1','SymLocal2', ...
                'SymLocal3','SymLocal4','SymLocal5','SymLocal6', ...
                'L1','L23Trans','L45Sym','L45Peak','L56Trans', ...
                'SymLocalaxon','SymLocal1axon','SymLocal2axon', ...
                'SymLocal3axon','SymLocal4axon','SymLocal5axon','SymLocal6axon',...
                'L1axon','L23Transaxon','L45Symaxon','L45Peakaxon','L56Transaxon'};
        case 'v2'
            CellTypeList = {'L2PY','L3PY','L4PY','L4sp','L4ss','L5IT', ...
                'L5PT','L6ACC','L6BCC','L6CT','INH','VPM'};
        otherwise
            CellTypeList = {'L2','L34','L4py','L4sp','L4ss','L5st', ...
                'L5tt','L6cc','L6ccinv','L6ct','SymLocal1','SymLocal2', ...
                'SymLocal3','SymLocal4','SymLocal5','SymLocal6','L1', ...
                'L23Trans','L45Sym','L45Peak','L56Trans'};
    end
end