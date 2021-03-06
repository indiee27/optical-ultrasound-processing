function alines = importfile(filename, dataLines)
%IMPORTFILE Import data from a text file
%  ALINES = IMPORTFILE(FILENAME) reads data from text file FILENAME for
%  the default selection.  Returns the data as a cell array.
%
%  ALINES = IMPORTFILE(FILE, DATALINES) reads data for the specified row
%  interval(s) of text file FILENAME. Specify DATALINES as a positive
%  scalar integer or a N-by-2 array of positive scalar integers for
%  dis-contiguous row intervals.
%
%  Example:
%  alines = importfile("C:\Users\where\Documents\Work\PhD things\Practice Data\M mode imaging data processing example\X1\alines.txt", [1, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 29-Sep-2020 12:43:23

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [1, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ":";

% Specify column names and types
opts.VariableNames = ["Date", "VarName2"];
opts.VariableTypes = ["char", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Date", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Date", "EmptyFieldRule", "auto");

% Import the data
alines = readtable(filename, opts);

%% Convert to output type
alines = table2cell(alines);
numIdx = cellfun(@(x) ~isnan(str2double(x)), alines);
alines(numIdx) = cellfun(@(x) {str2double(x)}, alines(numIdx));
end