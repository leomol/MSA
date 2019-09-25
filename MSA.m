%% Process miniscope data and quantify activity (e.g. spontaneous peaks and correlations) on arbitrary epoch definitions.
% 
% Epochs consists of labeled arrays of timestamps (s) and share a common window size.
% For example, 
%   epochs = {'A', 100, 'B', [200 2000], 'C', 300};
%   windowSize = 10;
% Means there are three epochs called A, B, and C:
%   A: 100 to 110
%   B: 200 to 210 and 2000 to 2010
%   C: 300 to 310

% 2019-03-22. Leonardo Molina.
% 2019-09-25. Last modified.

%% Add dependencies.
addpath('common');

%% Settings.
% epochs: Epochs contain start times relative to the beginning of the recording.
% window: Length of the windows to compare (s).
% dataFilename = 'C:/Users/molina/Documents/public/HALO/data/Miniscope/miniscope STS filtered.xlsx';
% data = loadData(dataFilename, 1);
% epochs = {'Pre', 60, 'During', 360, 'FS', 660, 'post' 960};
% epochWindow = 180;

dataFilename = 'C:\Users\molina\Documents\public\MATLAB\MSA\data\Miniscope.csv';
ttlFilename = 'C:\Users\molina\Documents\public\MATLAB\MSA\data\MiniscopeTTL.csv';
data = loadData(dataFilename);
ttl = loadDoricTTL(ttlFilename);
epochs = {'TTL stimulation', ttl, 'TTL baseline', ttl - 4};
epochWindow = 3;

% dataFilename = 'C:/Users/molina/Documents/public/HALO/data/P0001.xlsx';
% data = loadData(dataFilename, 2);
% epochs = {'Pre', 180, 'During', 420, 'Post', 720};
% epochWindow = 60;

% dataFilename = 'C:/Users/molina/Documents/public/HALO/data/Octavia/Octavia 17sec airpuff.csv';
% ttlFilename = 'C:/Users/molina/Documents/public/HALO/data/Octavia/Octavia 17sec airpuff TTL.csv';
% data = loadData(dataFilename);
% ttl = loadDoricTTL(ttlFilename);
% epochs = {'TTL stimulation', ttl, 'TTL baseline', ttl - 4};
% epochWindow = 3;

% Max time offset to compute cross-correlation (s).
lag = 0.5;
% Peak detection.
peakThreshold = 2.0;
peakThresholdingFunction = @mad;
peakLowpassFrequency = 5;
triggeredWindow = 2;
% Hierarchical clustering parameters.
clusteringDepth = 5;
clusteringPercentile = 5;
cmap = lines();

%% Derived terms.
time = data(:, 1);
data = data(:, 2:end);
nSamples = size(data, 1);
nCells = size(data, 2);
period = mean(diff(time));
frequency = 1 / period;
nEpochWindow = round(epochWindow * frequency);
maxLag = round(lag * frequency);
nLags = 2 * maxLag + 1;
lag0 = (nLags + 1) / 2;
nEpochs = numel(epochs) / 2;
eventCounts = cellfun(@numel, epochs(2:2:end));
eventStart = [0, cumsum(eventCounts(1:end - 1))] + 1;
eventEnd = eventStart + eventCounts - 1;

% Normalize traces.
data = zscore(data);

%% Find peaks.
filterOrder = 12;
peaksFilter = designfilt('lowpassiir', 'HalfPowerFrequency', peakLowpassFrequency, 'SampleRate', frequency, 'DesignMethod', 'butter', 'FilterOrder', filterOrder);
lowPassData = filtfilt(peaksFilter, data);
peakBool = false(size(data));
w = ceil(0.5 * frequency / peakLowpassFrequency);
if mod(w, 2) ~= 0
    w = w + 1;
end
warning('OFF', 'signal:findpeaks:largeMinPeakHeight');
for c = 1:nCells
    % Find peak in filtered data.
    threshold = mean(lowPassData(:, c)) + peakThreshold * peakThresholdingFunction(lowPassData(:, c));
    [~, pId] = findpeaks(lowPassData(:, c), 'MinPeakHeight', threshold);
    % Find peak around unfiltered data (within a time window compatible with the filter).
    for i = 1:numel(pId)
        range = max(1, pId(i) - w / 2):min(nSamples, pId(i) + w / 2);
        [~, k] = max(+data(range, c));
        pId(i) = k + range(1) - 1;
    end
    peakBool(pId, c) = true;
end

epochPeakCount = zeros(nEpochs, nCells);
epochPeaksBool = false(nSamples, nCells);
for e = 1:nEpochs
    timeLimits = epochs{2 * e};
    timeLimits(2, :) = timeLimits(1, :) + epochWindow;
    epochId = time2id(time, timeLimits);
    epochBool = false(nSamples, 1);
    epochBool(epochId) = true;
    partial = bsxfun(@and, epochBool, peakBool);
    epochPeaksBool = epochPeaksBool | partial;
    epochPeakCount(e, :) = sum(partial, 1);
end

%% Save peak counts per epoch and cell to file.
[folder, basename] = fileparts(dataFilename);
output = fullfile(folder, sprintf('%s - peak count.csv', basename));
fid = fopen(output, 'w');
fprintf(fid, 'Peak count. Columns = cells. Rows = epochs.\n');
fprintf(fid, [strjoin(repmat({'%d'}, 1, nCells), ',\t') '\n'], epochPeakCount');
fclose(fid);

%% Create time windows for each event.
t = cat(2, epochs{2:2:end});
nEvents = numel(t);
% Locate in time vector.
ids = ceil((t(:) + eps) / period);
% Indexed window for each id.
ids = bsxfun(@plus, (1:nEpochWindow) - 1, ids)';

% Remove out of range traces from triggered events.
bool = all(ids >= 1 & ids <= nSamples, 1);
ids = ids(:, bool);
if ~all(bool)
    fprintf(2, 'Removed %i event(s) because your window definition leads to out of range values.\n', sum(~bool));
end

% Use offset to replicate indices for every cell.
offset = ones(nEpochWindow, nEvents, nCells);
offset(:, :, 1) = 0;
offset = nSamples * cumsum(offset, 3);
ids = bsxfun(@plus, ids, offset);

%% Pair-wise cross-correlation.
ticker = tic;
pairs = nchoosek(1:nCells, 2);
nPairs = size(pairs, 1);
xcAll = NaN(nLags, nEvents, nPairs);
for e = 1:nEvents
    % Slice for multi-threading.
    traces1 = zeros(nEpochWindow, nPairs);
    traces2 = zeros(nEpochWindow, nPairs);
    for p = 1:nPairs
        p1 = pairs(p, 1);
        p2 = pairs(p, 2);
        traces1(:, p) = data(ids(:, e, p1));
        traces2(:, p) = data(ids(:, e, p2));
    end
    parfor p = 1:nPairs
        xcAll(:, e, p) = xcorr(traces1(:, p), traces2(:, p), maxLag);
    end
end
% Normalize.
xcAll(:) = xcAll(:) / nEpochWindow;
fprintf('%d pair-wise cross-correlations (%d cells x %d events x %d lags) in %.2fs.\n', nPairs * nEvents * nLags, nCells, nEvents, nLags, toc(ticker));

%% Plot traces and event windows.
name = 'Traces and events';
figure('name', name);
% Add space between traces.
pad = 5 * std(data(:));
separation = pad * ones(1, nCells);
separation(1) = 0;
separation = cumsum(separation);
ylims = [separation(1) + min(data(:, 1)), separation(end) + max(data(:, end))];
patches = cell(nEpochs, 1);
for e = 1:nEpochs
    epochTimes = epochs{2 * e};
    epochName = epochs{2 * e - 1};
    [faces, vertices] = patchEpochs([epochTimes; epochTimes + epochWindow], ylims(1), ylims(2));
    patches{e} = patch('Faces', faces, 'Vertices', vertices, 'FaceColor', cmap(e, :), 'EdgeColor', 'none', 'FaceAlpha', 0.50, 'DisplayName', sprintf('%s | peaks:%i', epochName, sum(epochPeakCount(e, :))));
end
patches = [patches{:}];
hold('all');
xx = repmat(time, 1, nCells);
yy = bsxfun(@plus, separation, data);
hs = plot(xx, yy, 'LineWidth', 0.5, 'HandleVisibility', 'off');
plot(xx(epochPeaksBool), yy(epochPeaksBool), 'k.', 'HandleVisibility', 'off');
axis('tight');
legend('show');
set(gca, 'YTick', separation, 'YTickLabel', 1:nCells);
xlabel('Time (s)');

%% Verify ids by plotting one trace at a time.
h = plot(NaN(2, 1), NaN(2, 1), 'LineWidth', 2, 'LineStyle', '-', 'Color', [1, 0, 0], 'HandleVisibility', 'off');
for e = 1:nEvents
    for c = 1:nCells
        set(h, 'XData', time(ids(:, e, 1)), 'YData', data(ids(:, e, c)) + separation(c));
        pause(0.025);
    end
end
delete(h);

%% Verify cross-correlation for a given pair of cells.
pair = [2, 10];
pairId = find(all(pairs == pair(1) | pairs == pair(2), 2));
set(hs, 'LineWidth', 0.5);
set(hs(pair), 'LineWidth', 4);
for e = 1:nEpochs
    epochName = epochs{2 * e - 1};
    xc = xcorr(data(ids(:, e, pair(1))), data(ids(:, e, pair(2))), 0) / nEpochWindow;
    patches(e).DisplayName = sprintf('xcorr:%.4f==%.4f? | %s | peaks:%i.', xcAll(lag0, e, pairId), xc, epochName, sum(epochPeakCount(e, :)));
end

%% Plot cross-correlation.
name = 'Cross-correlation';
figure('name', name);
hold('all');
t = linspace(-0.5 * lag, 0.5 * lag, nLags);
for e = 1:nEpochs
    epochName = epochs{2 * e - 1};
    xc = mean(xcAll(:, eventStart(e):eventEnd(e), :), 3);
    av = mean(xc, 2)';
    se = std(xc, [], 2)' / sqrt(eventCounts(e));
    plot(t, av, 'Color', cmap(e, :), 'LineStyle', '-', 'DisplayName', epochName);
    faces = 1:2 * nLags;
    vertices = [t; av + se / 2];
    vertices = cat(2, vertices, [fliplr(t); fliplr(av - se / 2)])';
    patch('Faces', faces, 'Vertices', vertices, 'FaceColor', cmap(e, :), 'EdgeColor', 'none', 'FaceAlpha', 0.10, 'HandleVisibility', 'off');
end
title(name);
xlabel('time (s)');
ylabel('Cross-correlation');
legend('show');

%% Plot cross-correlation (zero lag) for each epoch.
name = 'Cross-correlation (zero lag)';
figure('name', name);
title(name);
hold('all');
for e = 1:nEpochs
    epochTimes = epochs{2 * e};
    epochName = epochs{2 * e - 1};
    xc = mean(mean(xcAll(lag0, eventStart(e):eventEnd(e), :), 3), 2);
    plot(epochTimes, xc, 'Color', cmap(e, :), 'LineStyle', '-', 'Marker', 'o', 'DisplayName', epochName);
end
xlabel('Time (s)');
ylabel('Cross-correlation');
legend('show');

%% Plot cross-correlation matrix at zero lag.
name = 'Cross-correlation matrix at zero lag. Pairs sorted to satisfy hierarchical clustering.';
figure('name', name);
xc = xcAll(lag0, :, :);
clims = [min(xc(:)), max(xc(:))];
for e = 1:nEpochs
    epochName = epochs{2 * e - 1};
    xc = squeeze(mean(xcAll(lag0, eventStart(e):eventEnd(e), :), 2));
    xcMatrix = squareform(xc);
    xcMatrix(eye(nCells, 'logical')) = clims(2);
    if e == 1
        z = linkage(xc(:)', 'weighted');
        coef = inconsistent(z, clusteringDepth);
        cutoff = prctile(coef(coef(:, 4) > 0, 4), clusteringPercentile);
        clusters = cluster(z, 'criterion', 'distance', 'cutoff', cutoff);
        [~, order] = sort(clusters);
    end
    subplot(1, nEpochs, e);
    imagesc(xcMatrix(order, order));
    axis('equal', 'square', 'tight')
    caxis(clims);
    title(epochName);
end