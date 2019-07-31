function [data, names, sheetName] = loadInscopixData(filename, number)
    [~, sheetNames] = xlsfinfo(filename);
    if nargin == 1
        number = 1;
    end
    sheetName = sheetNames{number};
    state = warning('QUERY', 'MATLAB:table:ModifiedAndSavedVarnames');
    warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');
    table = readtable(filename, 'Sheet', sheetName);
    warning(state.state, 'MATLAB:table:ModifiedAndSavedVarnames');
    data = table2array(table);
    names = table.Properties.VariableNames;
end