% [id, limits] = time2id(time, epochs)
% Returns 
% id: index of each time value enclosed by epoch pairs (from1, to1, from2, to2, ...).
% limits: index of each epoch pair.

% 2019-02-01. Leonardo Molina.
% 2022-06-20. Last modified.
function [ids, limits] = time2id(time, epochs)
    epochs = epochs(:);
    % timeLimits = zeros(2, nEpochs);
    a = arrayfun(@(t) find(time >= t, 1, 'first'), epochs(1:2:end), 'UniformOutput', false);
    b = arrayfun(@(t) find(time <= t, 1, 'last'), epochs(2:2:end), 'UniformOutput', false);
    k = cellfun(@isempty, a) | cellfun(@isempty, b);
    a(k) = [];
    b(k) = [];
    nEpochs = numel(a);
    limits = zeros(2, nEpochs);
    limits(1:2:end) = [a{:}];
    limits(2:2:end) = [b{:}];
    % When "last" can't find anything ahead of "first", force a single time point.
    k = diff(limits, [], 1) < 0;
    limits(2, k) = limits(1, k);
    ids = arrayfun(@(e) colon(limits(1, e), limits(2, e))', 1:nEpochs, 'UniformOutput', false);
    ids = cat(1, ids{:});
end