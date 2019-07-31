function id = time2id(time_vector, time_limits)
    time_limits = time_limits(:);
    nEpochs = numel(time_limits) / 2;
    epochs = zeros(2, nEpochs);
    epochs(1:2:end) = arrayfun(@(l) find(time_vector >= l, 1, 'first'), time_limits(1:2:end));
    epochs(2:2:end) = arrayfun(@(l) find(time_vector <= l, 1, 'last'), time_limits(2:2:end));
    id = arrayfun(@(e) epochs(1, e):epochs(2, e), 1:nEpochs, 'UniformOutput', false);
    id = [id{:}];
end