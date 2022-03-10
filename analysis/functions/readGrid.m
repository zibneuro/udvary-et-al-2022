function tbl = readGrid(filename)
%% Reads subvolume grid and returns as table
    dataLines = [2, Inf];
    opts = delimitedTextImportOptions("NumVariables", 6);
    opts.DataLines = dataLines;
    opts.Delimiter = ",";
    opts.VariableNames = ["ix", "iy", "iz", ...
                    "cube_origin_x", "cube_origin_y", "cube_origin_z"];
    opts.VariableTypes = ["double", "double", "double", ...
                    "double", "double", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    tbl = readtable(filename, opts);

end