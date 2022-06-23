% [data, names] = loadData(filename)
% Load one of several formats of data expected from a recording with either:
% Doric (csv), Multifiber (csv), Inscopix (csv), XLS files, or Axon (abf).
% 
% [data, names] = loadData(filename, sheetName)
% [data, names, sheetName] = loadData(filename, sheetNumber)
% If the file has multiple sheets, use sheetName or sheetNumber to select one.

% 2019-05-07. Leonardo Molina.
% 2022-06-23. Last modified.
function [data, names, sheetName] = loadData(filename, varargin)
    [~, ~, extension] = fileparts(filename);
    if ismember(lower(extension), {'.xls', '.xlsx'})
        [data, names, sheetName] = loadXLS(filename, varargin{:});
    elseif ismember(lower(extension), '.abf')
        [data, ~, names] = loadABF(filename, varargin{:});
        sheetName = '';
    else
        [data, names] = loadCSV(filename);
        sheetName = '';
    end
end

function [data, names, sheetName] = loadXLS(filename, sheet)
    if nargin == 1
        sheet = 1;
    end
    if isnumeric(sheet)
        [~, sheetNames] = xlsfinfo(filename);
        sheetName = sheetNames{sheet};
    else
        sheetName = sheet;
    end
    state = warning('QUERY', 'MATLAB:table:ModifiedAndSavedVarnames');
    warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');
    table = readtable(filename, 'Sheet', sheetName);
    warning(state.state, 'MATLAB:table:ModifiedAndSavedVarnames');
    data = table2array(table);
    names = table.Properties.VariableNames;
end

function [data, names] = loadCSV(filename)
    % Read first few lines to determine type case.
    fid = fopen(filename, 'r');
    lines = cell(3, 1);
    for l = 1:numel(lines)
        line = fgetl(fid);
        line = textscan(line, '%s', 'Delimiter', ',');
	    lines{l} = line{1}(:)';
    end

    is1Number = cellfun(@validateNumber, lines{1});
    is2Number = cellfun(@validateNumber, lines{2});
    is3Number = cellfun(@validateNumber, lines{3});
    is1Name = ~(is1Number | cellfun(@isempty, lines{1}));
    is2Name = ~(is2Number | cellfun(@isempty, lines{2}));
    firstIsBlank = isempty(strtrim(lines{1}{1}));
    hasKeepLine = is2Name(1) && all(cellfun(@(x) ismember(x, {'1', '0', 'true', 'false'}), lower(lines{2}(2:end))));
    nColumns = numel(lines{3});

    success = true;
    if all(is1Number) && all(is2Number) && all(is3Number)
        % "Standard" variation 1
        %   number0, number1, ...
        %   number0, number1, ...
        %   number0, number1, ...
        keep = true(1, nColumns);
        names = ['time', arrayfun(@(x) sprintf('c%i', x), 1:nColumns, 'UniformOutput', false)];
        nHeaderLines = 0;
    elseif all(is1Name) && all(is2Number)
        % "Standard" variation 2
        %   header0, header1, header2, ...
        %   number0, number1, number2, ...
        keep = true(1, nColumns);
        names = lines{1};
        nHeaderLines = 1;
    elseif firstIsBlank && all(is1Name(2:end)) && all(is2Number) && all(is3Number)
        % "Standard" variation 3
        %   blank,   header1, header2, ...
        %   number0, number1, number2, ...
        %   number0, number1, number2, ...
        keep = true(1, nColumns);
        names = ['time', lines{1}(2:end)];
        nHeaderLines = 1;
    elseif is1Name(1) && all(cellfun(@isempty, lines{1}(2:end))) && all(is2Number) && all(is3Number)
        % "Segmentation"
        %   header0 <,blank>...
        %   number0, number1, ...
        %   number0, number1, ...
        keep = true(1, nColumns);
        names = [lines{1}{1}, arrayfun(@(x) sprintf('c%i', x), 1:nColumns - 1, 'UniformOutput', false)];
        nHeaderLines = 1;
    elseif firstIsBlank && all(is1Name(2:end)) && all(is3Number)
        if hasKeepLine
            % "nVista" variation 1
            %   blank,   header1, header2, ...
            %   header0,   keep1,   keep2, ...
            %   number0, number1, number2, ...
            keep = [true, cellfun(@(x) ismember(x, {'1', 'true'}), lower(lines{2}(2:end)))];
            names = ['time', lines{1}(2:end)];
            nHeaderLines = 2;
        elseif all(is2Name)
            % "nVista" variation 2
            %   blank,   header1, header2, ...
            %   header0, header1, header2, ...
            %   number0, number1, number2, ...
            keep = true(1, nColumns);
            names = ['time', lines{1}(2:end)];
            nHeaderLines = 2;
        else
            success = false;
        end
    elseif all(is1Name) && all(is2Number == is3Number)
        % "Multifiber"
        %   header0, header1, header2, ...
        %   number0, number1, ... <not-number>, <number> ...
        %   number0, number1, ... <not-number>, <number> ...
        keep = is2Number;
        names = lines{1};
        nHeaderLines = 1;
    else
        success = false;
    end

    if success
        % Reset carret.
        fseek(fid, 0, 'bof');
        code = 'sf';
        format = [repmat('%', 1, nColumns); code(is3Number + 1)];
        data = textscan(fid, format(:)', 'Delimiter', ',', 'HeaderLines', nHeaderLines);
        fclose(fid);
        data = cat(2, data{keep});
        names = names(keep);
        % Remove columns with only NaN values.
        c = all(isnan(data), 1);
        data(:, c) = [];
        names(c) = [];
        % Remove rows with any NaN values.
        r = any(isnan(data), 2);
        data(r, :) = [];
    else
        error('Error loading data: Unknown format');
    end
end

function result = validateNumber(numbers)
    result = ~isnan(str2double(numbers)) | contains(numbers, 'nan', 'IgnoreCase', true);
end