function tbl = readCellTypesCSV(filename)
%% Reads CellTypes csv file and returns as table
    dataLines = [2, Inf];
    opts = delimitedTextImportOptions("NumVariables", 3);
    opts.DataLines = dataLines;
    opts.Delimiter = ",";
    opts.VariableNames = ["celltypeID", "name", "excitatory"];
    opts.VariableTypes = ["double", "string", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts = setvaropts(opts, "name", "WhitespaceRule", "preserve");
    opts = setvaropts(opts, "name", "EmptyFieldRule", "auto");
    tbl = readtable(filename, opts);
end