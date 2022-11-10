function TG100119720210601T120000 = ReadLI7810(filename, dataLines)
%IMPORTFILE Import data from a text file
%  TG100119720210601T120000 = IMPORTFILE(FILENAME) reads data from text
%  file FILENAME for the default selection.  Returns the data as a table.
%
%  TG100119720210601T120000 = IMPORTFILE(FILE, DATALINES) reads data for
%  the specified row interval(s) of text file FILENAME. Specify
%  DATALINES as a positive scalar integer or a N-by-2 array of positive
%  scalar integers for dis-contiguous row intervals.
%
%  Example:
%  TG100119720210601T120000 = importfile("C:\Users\thmorin\OneDrive - SUNY ESF\Documents\Projects\WillowGHG\data\TG10-01197-2021-06-01T120000.data", [8, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 08-Jun-2021 10:16:50

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [8, Inf];
end

%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 21);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["DATAH", "SECONDS", "NANOSECONDS", "NDX", "DIAG", "DATE", "TIME", "H2O", "CO2", "CH4", "CAVITY_P", "CAVITY_T", "LASER_PHASE_P", "LASER_T", "RESIDUAL", "RING_DOWN_TIME", "THERMAL_ENCLOSURE_T", "PHASE_ERROR", "LASER_T_SHIFT", "INPUT_VOLTAGE", "CHK"];
opts.VariableTypes = ["categorical", "double", "double", "double", "double", "datetime", "datetime", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts = setvaropts(opts, 6, "InputFormat", "yyyy-MM-dd");
opts = setvaropts(opts, 7, "InputFormat", "HH:mm:ss");
opts = setvaropts(opts, 1, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
TG100119720210601T120000 = readtable(filename, opts);

end