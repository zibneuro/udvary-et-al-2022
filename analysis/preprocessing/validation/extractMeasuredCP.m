function [] = extractMeasuredCP(pathToMeasuredExcelTable,pathToStoreMatFiles)
% extractMeasuredCP(pathToMeasuredExcelTable,pathToStoreMatFiles)
% Save experimental measurements in pathToStoreMatFiles/dataComparisons\MeasuredCP.mat
% Input:
% - pathToMeasuredExcelTable: location of .xlsx table
% - pathToStoreMatFiles: folder/to/save MeasuredCP.mat
% Output:
% - stores file containing
%     Authors [First Author et al.]
%     Year [year of publication]
%     Title [First words of Title]
%     PresynapseLayer [Layer or Presynaptic Cell]
%     PresynapseType [Presynaptic Cell Type]
%     PostsynapseLayer [Layer or Ppstsynaptic Cell]
%     PostsynapseType [Postsynaptic Cell Type]
%     CP_min [minimum of measured Connection probability]
%     CP_max [maximum of measured Connection probability]
%     NumSamples [tested number of pairs]
%     SpeciesArea [Species and Area measured was conducted]
%     SliceWidth [Width of Slice in um]
%     SliceOrientation [Orientation of Slice]
%     ISLimit_y [Range Intersomatic Distance in y] [1 x 2]
%     ISLimit_z [Range Intersomatic Distance in z] [1 x 2]
%     ISLimit_yz [Range Intersomatic Distance in yz] [1 x 2]
%     ISLimit_xyz [Range Intersomatic Distance in xyz] [1 x 2]
%     TissueDepthLimit [Range of tissue depth of cell pair] [1 x 2]
%     InsideColumn [Cell pairs insideColumn (1), in septum (0), or both
%       [nan]
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)

    matlabPath = 'D:\udvary-et-al-2022\analysis\';

    if nargin<1
        pathToMeasuredExcelTable = [matlabPath 'preprocessing\data\MeasuredCP.xlsx'];
    end
    if nargin<2
        pathToStoreMatFiles = [matlabPath 'data\'];
    end
    
    array = xlsread(pathToMeasuredExcelTable);
    [~, ~, raw] = xlsread(pathToMeasuredExcelTable,'Sheet1', ...
                                ['A2:AD' num2str(size(array,1)+1)]);
    raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
    cellVectors = raw(:,[1,3:7,11:14,26:30]); % Strings
    raw = raw(:,[2,8:10,15:25]); % Numbers
    raw(cellfun(@(x)~isnumeric(x),raw)) = {nan};
    data = reshape([raw{:}],size(raw));

    Authors = cellVectors(:,1);
    Title = cellVectors(:,2);
    PresynapseLayer = cellVectors(:,3);
    PresynapseType = cellVectors(:,4);
    PostsynapseLayer = cellVectors(:,5);
    PostsynapseType = cellVectors(:,6);
    Species = cellVectors(:,7);
    Area = cellVectors(:,8);
    SliceWidth = cellVectors(:,9);
    SliceOrientation = cellVectors(:,10);
    Synapses_Anat = strcmp(cellVectors(:,11),'yes');
    Synapses_Funct = strcmp(cellVectors(:,12),'yes');
    SharpElectrode = strcmp(cellVectors(:,13),'yes');
    ExpmtlMethod = cellVectors(:,14);
    AnimalAge = cellVectors(:,15);
    
    Year = data(:,1);
    CP_min = data(:,2);
    CP_max = data(:,3);
    NumSamples = data(:,4);
    ISLimit_y = [data(:,5) data(:,6)];
    ISLimit_z = [data(:,7) data(:,8)];
    ISLimit_yz = [data(:,9) data(:,10)];
    ISLimit_xyz = [data(:,11) data(:,12)];
    TissueDepthLimit = [data(:,13) data(:,14)];
    InsideColumn = data(:,15);    

    for i = 1:length(SliceWidth)
       if isnumeric(SliceWidth{i}) 
          SliceWidth{i} = num2str(SliceWidth{i});
       end
    end

    % Clear temporary variables
    clearvars array data raw cellVectors R i;
    save([pathToStoreMatFiles 'MeasuredCP.mat']);

end