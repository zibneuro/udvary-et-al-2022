function tbl = readNeuronsCSV(filename)
%% Reads Neurons.csv file (Neuron meta file) and returns as table
    dataLines = [2, Inf];

    %% Set up the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 12);

    % Specify range and delimiter
    opts.DataLines = dataLines;
    opts.Delimiter = ",";

    % Specify column names and types
    opts.VariableNames = ["id", "graph_id", "soma_x", "soma_y", "soma_z", ...
                        "cell_type", "nearest_column", "region", ...
                        "laminar_location", "cortical_depth", ...
                        "synaptic_side", "inside_vS1"];
    opts.VariableTypes = ["double", "double", "double", "double", ...
                        "double", "double", "double", "double", "double", ...
                        "double", "double", "double"];

    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Import the data
    tbl = readtable(filename, opts);

end