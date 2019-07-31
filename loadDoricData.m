function [data, names] = loadDoricData(filename)
    fid = fopen(filename, 'r');
    names = fgetl(fid);
    names = textscan(names, '%s', 'Delimiter', ',');
    names = names{1}(2:end);
    keep = fgetl(fid);
    keep = textscan(keep, '%s', 'Delimiter', ',');
    keep = ismember(keep{1}(2:end), 'accepted');
    data = textscan(fid, repmat('%f', 1, numel(names) + 1), 'Delimiter', ',', 'HeaderLines', 2);
    data = cat(2, data{:});
    names = names(keep);
    fclose(fid);
end