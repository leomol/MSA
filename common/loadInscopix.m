% [data, cellIds] = loadInscopix(filename)
% Return data from an Inscopix csv file with only columns marked as 'accepted'.

% 2023-07-19. Leonardo Molina.
% 2023-07-19. Last modified.
function [data, cellIds] = loadInscopix(filename)
    fid = fopen(filename, 'r');
    line = fgetl(fid);
    columns = regexp(line, 'C(\d+)', 'tokens');
    cellIds = str2double([columns{:}]);
    line = fgetl(fid);
    items = strsplit(line, ', ');
    items = items(2:end);
    accepted = ismember(items, 'accepted');
    nColumns = numel(accepted) + 1;
    
    nHeaderLines = 2;
    fseek(fid, 0, 'bof');
    format = repmat('%f', 1, nColumns);
    data = textscan(fid, format(:)', 'Delimiter', ',', 'HeaderLines', nHeaderLines);
    fclose(fid);
    data = cat(2, data{[true, accepted]});
end