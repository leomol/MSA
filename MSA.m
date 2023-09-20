% msa = MSA(data, configuration);
% 
% Normalize, filter, and detect peaks of spontaneous
% activity in user defined epochs.
% 
% Each column in data corresponds to activity from individual cells.
% 
% configuration is a struct with the following fields (defaults are used for missing fields):
%     conditionEpochs - Epochs for different conditions: {'epoch1', [start1, end1, start2, end2, ...], 'epoch2', ...}
%     thresholdEpochs - Epochs to include for peak threshold calculation.
%     events - Event-triggered data; times at which a type of event occurs.
%     resamplingFrequency - Resampling frequency (Hz).
%     lowpassFrequency - Lowest frequency permitted in normalized signal.
%     peakSeparation - Minimum time separation between two peaks.
% 
% Processing steps:
%   -Resample signal and reference to a given frequency.
% 	-Baseline correction modeled as an exponential decay of the low-pass filtered data (optionally using airPLS).
% 	-Correct for motion artifacts by subtracting reference to signal, after a polynomial fit (optional).
%   -Remove fast oscillations with a low-pass filter.
% 	-Normalize data as df/f or z-score according to settings.
% 	-Find peaks of spontaneous activity in low-pass filtered data.
% 
% Normalization is calculated as (f - f0) / f1 where f0 and f1 can be data provided by the
% user or calculated using given functions:
% 
%     Normalization from given functions:
%         f0 and f1 are common to all data points and calculated from all data:
%             df/f:
%                 configuration.f0 = @mean;
%                 configuration.f1 = @mean;
%             z-score:
%                 configuration.f0 = @mean;
%                 configuration.f1 = @std;
%             z-score - option 1:
%                 configuration.f0 = @median;
%                 configuration.f1 = @mad;
%             z-score - option 2:
%                 configuration.f0 = @median;
%                 configuration.f1 = @std;
% 
%         f0 and f1 are common to all data points and calculated at given epochs:
%             df/f:
%                 epochs = [0, 100, 500, 550, 1000, Inf]
%                 configuration.f0 = {@mean, epochs};
%                 configuration.f1 = {@mean, epochs};
% 
%         f0 and f1 are calculated for each data point based on a moving window:
%             df/f:
%                 window = 60;
%                 configuration.f0 = {@movmean, window};
%                 configuration.f1 = {@movmean, window};
%             (further combinations possible with @movmean, @movmedian, @movstd, @movmad, @mov...).
% 
%     Normalization from given data:
%         f0 = ones(size(time));
%         f1 = ones(size(time)) * 10;
%         configuration.f0 = f0;
%         configuration.f1 = f1;
% 
% Fluorescence deflections are considered peaks when they exceed a threshold calculated as
% k * f2 + f3 and they are provided by the user as configuration.threshold = {k, f2, f3}
% Examples:
%   2.91 median absolute deviations from the median:
%     configuration.threshold = {2.91, @mad, @median}
%   2.91 median absolute deviations from 0:
%     configuration.threshold = {2.91, @mad, 0}
%   2.00 standard deviations from the mean:
%     configuration.threshold = {2.00, @std, @mean}
% 
% See examples and source code for detailed analysis steps and default parameters.
% Units for time and frequency are seconds and hertz respectively.
% 
% 2019-03-22. Leonardo Molina.
% 2023-09-20. Last modified.
classdef MSA < handle
    properties
        configuration
        
        warnings
        time
        frequency
        peakIds
        peakLabels
        eventIds
        eventLabels
        epochIds
        epochLabels
        peakCounts
        spikeCounts
        
        signalRaw
        signalBaseline
        signal
        f
        fSmooth
        f0
        f1
        dff
        spikes
        duration
        area
        peakThreshold

        muaTime
        muaCounts
    end
    
    properties (Access = private)
        nConditions
        nCells
        cleanId
        windowTemplate
        % Peaks detected including artifacts.
        peakWidth
        peakMaskAll
        normalizedArea
        firingRate
        boxplotIds
        boxplotLabels
        spikePeakIds

        % Settings for visualization.
        cmap = lines();
        zoomSettings = {0.99, 0.50};
        signalColor = [0, 0.4470, 0.7410];
        alternateColor = [0, 0.6470, 0.9410];
        peaksLineColor = [0.4660, 0.6740, 0.1880];
        peaksMarkerColor = [1, 0, 0];
        dashColor = [0, 0, 0];
    end

    methods
        function obj = MSA(data, parameters)
            % Time is in the first column and is common to all data.
            time = data(:, 1);
            % Separate signal from time.
            signal = data(:, 2:end);
            
            % Accumulate all warnings.
            obj.warnings = {};
            
            % The variable parameters is optional.
            if nargin < 2
                parameters = struct();
            end
            
            % Use defaults for missing parameters.
            defaults.artifactEpochs = [];
            defaults.conditionEpochs = {'Data', [-Inf, Inf]};
            defaults.events = [];
            defaults.resamplingFrequency = NaN;
            defaults.lowpassFrequency = Inf;
            defaults.peakSeparation = 0.5;
            defaults.f0 = 0.0;
            defaults.f1 = 1.0;
            defaults.thresholdEpochs = NaN;
            defaults.threshold = {2.91, @mad, @median};
            defaults.triggeredWindow = 10.0;
            defaults.spikes = [];
            defaults.peakDetectionMode = 'height';
            
            % Override defaults with user parameters.
            configuration = defaults;
            defaultNames = fieldnames(defaults);
            parametersNames = fieldnames(parameters);
            for i = 1:numel(parametersNames)
                name = parametersNames{i};
                if ismember(name, defaultNames)
                    configuration.(name) = parameters.(name);
                else
                    obj.warnings{end + 1} = warn('[parsing] "%s" is not a valid parameter.', name);
                end
            end
            
            options = {'MinPeakHeight', 'MinPeakProminence'};
            k = strcmpi(configuration.peakDetectionMode, {'height', 'prominence'});
            if any(k)
                peakDetectionMode = options{k};
            else
                error('peakDetectionMode must be either "height" or "prominence"');
            end
            
            % Resampling frequency defaults to the smallest between 100Hz and the source frequency.
            sourceFrequency = 1 / median(diff(time));
            if isnan(configuration.resamplingFrequency)
                configuration.resamplingFrequency = min(100, sourceFrequency);
            end
            
            % Resample to target frequency.
            if configuration.resamplingFrequency < sourceFrequency
                frequency = configuration.resamplingFrequency;
                % Express frequency as a ratio p/q.
                [p, q] = rat(frequency / sourceFrequency);
                % Resample: interpolate every p/q/f, upsample by p, filter, downsample by q.
                [signal, time2] = resample(signal, time, frequency, p, q);
                time = time2;
            elseif configuration.resamplingFrequency > sourceFrequency
                frequency = sourceFrequency;
                obj.warnings{end + 1} = warn('[resampling] Cannot resample to frequencies higher than the source frequency (%.2f Hz).', sourceFrequency);
            else
                frequency = sourceFrequency;
            end
            
            % Setup.
            [nSamples, nCells] = size(signal);
            nConditions = numel(configuration.conditionEpochs) / 2;
            
            % Index of artifacts and non-artifacts.
            allIds = colon(1, nSamples)';
            artifactId = time2id(time, configuration.artifactEpochs);
            artifactFreeId = setdiff(allIds, artifactId);
            
            % Define clean epochs for fitting and peak detection.
            % Remove data from the edges as filtering will introduce artifacts here.
            % The exclusion window is the smallest between half a second and 5% of sample size.
            edgeWindow =  min(ceil(0.5 * frequency / configuration.lowpassFrequency), round(0.05 * nSamples));
            edgeId = union(artifactId, [1:edgeWindow, numel(time) - edgeWindow + 1:nSamples]');
            edgeFreeId = setdiff(allIds, edgeId);
            cleanId = intersect(artifactFreeId, edgeFreeId);
            
            % Threshold epochs defaults to everything.
            if isnan(configuration.thresholdEpochs)
                thresholdId = cleanId;
            else
                thresholdId = time2id(time, configuration.thresholdEpochs);
                thresholdId = intersect(cleanId, thresholdId);
            end
            
            % Low-pass filter.
            if configuration.lowpassFrequency == Inf
                fSmooth = signal;
            else
                fFilter = designfilt('lowpassiir', 'HalfPowerFrequency', configuration.lowpassFrequency, 'SampleRate', frequency, 'DesignMethod', 'butter', 'FilterOrder', 12);
                fSmooth = filtfilt(fFilter, signal);
            end
            
            % Normalize.
            f0 = parseNormalization(configuration.f0, fSmooth, time);
            f1 = parseNormalization(configuration.f1, fSmooth, time);
            dff = (fSmooth - f0) ./ f1;
            
            % Get peak threshold.
            state = warning('Query', 'signal:findpeaks:largeMinPeakHeight');
            warning('Off', 'signal:findpeaks:largeMinPeakHeight');
            
            peakWidth = zeros(0, 1);
            peakMaskAll = false(size(dff));
            peakThreshold = NaN(nCells, 1);
            % Ids for all conditions, for a single cell.
            ids = time2id(time, [configuration.conditionEpochs{2:2:end}]);
            for u = 1:nCells
                peakThreshold(u) = threshold(configuration.threshold, dff(thresholdId, u));
                if strcmp(peakDetectionMode, 'MinPeakProminence')
                    peakThreshold(u) = abs(peakThreshold(u));
                end
                if any(+dff(:, u) >= peakThreshold(u))
                    [~, k, peakWidthCell] = findpeaks(+dff(:, u), peakDetectionMode, peakThreshold(u), 'MinPeakDistance', configuration.peakSeparation * frequency, 'WidthReference', 'halfheight');
                    peakWidthCell = peakWidthCell / frequency;
                    peakMaskAll(k, u) = true;
                    peakWidth = cat(1, peakWidth, peakWidthCell(ismember(k, ids)));
                end
            end
            
            warning(state.state, 'signal:findpeaks:largeMinPeakHeight');
            [~, halfWindow] = forceOdd(configuration.triggeredWindow * frequency);
            % Index template to apply around each peak.
            obj.windowTemplate = -halfWindow:halfWindow;
            
            % Get spiking statistics.
            spikePeakIds = false(nSamples, nCells);
            if isa(configuration.spikes, 'function_handle')
                fcn = configuration.spikes;
                spikes = zeros(nSamples, nCells);
                for u = 1:nCells
                    spikes(:, u) = fcn(dff(:, u));
                    [~, k] = findpeaks(spikes(:, u));
                    spikePeakIds(k, u) = true;
                end
            else
                spikes = configuration.spikes;
                configuration = rmfield(configuration, 'spikes');
            end

            % Get indices for epochs.
            % Start and stop vector indices for all provided epochs.
            epochIds = zeros(2, 0);
            % Numeric label corresponding to each epoch range.
            epochLabels = zeros(0, 1);
            % Vector index for each peak / event.
            peakIds = zeros(0, 1);
            eventIds = zeros(0, 1);
            % Numeric label corresponding to each peak / event.
            peakLabels = zeros(0, 1);
            eventLabels = zeros(0, 1);
            % Misc indexing / labeling.
            boxplotIds = zeros(0, 1);
            boxplotLabels = zeros(0, 1);
            peakCounts = zeros(nConditions, nCells);
            duration = zeros(nConditions, 1);
            area = zeros(nConditions, nCells);
            spikeCounts = zeros(nConditions, nCells);
            
            % Get indices for time triggers.
            x = arrayfun(@(t) find(time >= t, 1, 'first'), configuration.events, 'UniformOutput', false);
            k = ~cellfun(@isempty, x);
            eventTimeIds = [x{k}];
            
            linearIds = zeros(nSamples, nCells);
            linearIds(:) = 1:nSamples * nCells;
            
            for c = 1:nConditions
                % Accumulate vector indices limited to conditions.
                [ids, bounds] = time2id(time, configuration.conditionEpochs{2 * c});
                % Mask for this condition, for all cells.
                mask = repmat(ismember(1:nSamples, ids)', 1, nCells);

                % Start/stop-triggered data.
                n = numel(bounds) / 2 * nCells;
                epochIds = cat(2, epochIds, bounds);
                epochLabels = cat(1, epochLabels, repmat(c, n, 1));
                
                % Event-triggered data.
                eventIdsEpoch = linearIds(intersect(eventTimeIds, ids), :);
                n = numel(eventIdsEpoch);
                eventIds = cat(1, eventIds, eventIdsEpoch(:));
                eventLabels = cat(1, eventLabels, repmat(c, n, 1));
                
                % Peak-triggered data.
                peakIdsEpoch = find(mask & peakMaskAll);
                n = numel(peakIdsEpoch);
                peakIds = cat(1, peakIds, peakIdsEpoch);
                peakLabels = cat(1, peakLabels, repmat(c, n, 1));
                
                peakCounts(c, :) = sum(peakMaskAll(ids, :));
                if ~isempty(spikes)
                    spikeCounts(c, :) = sum(spikePeakIds(ids, :));
                end
                boxplotIds = cat(1, boxplotIds, ids);
                boxplotLabels = cat(1, boxplotLabels, repmat(c, numel(ids), 1));
                
                duration(c) = numel(ids) / frequency;
                area(c, :) = trapz(dff(ids, :)) / frequency;
            end
            
            % Replicate epoch ids for each cell column.
            epochIds = arrayfun(@(k) epochIds + k, linearIds(1, :) - 1, 'UniformOutput', false);
            epochIds = cat(1, epochIds{:});
            
            % Normalize area according to epoch length.
            normalizedArea = bsxfun(@rdivide, area, duration);
            normalizedArea(duration == 0) = 0;
            
            % Compute firing rate.
            firingRate = bsxfun(@rdivide, spikeCounts, duration);
            firingRate(duration == 0) = 0;
            
            obj.configuration = configuration;
            obj.time = time;
            obj.frequency = frequency;
            
            % Order depends on epoch definitions. Overlapping is possible and allowed.
            obj.peakWidth = peakWidth;
            obj.peakIds = peakIds;
            obj.peakLabels = peakLabels;
            obj.eventIds = eventIds;
            obj.eventLabels = eventLabels;
            obj.epochIds = epochIds;
            obj.epochLabels = epochLabels;
            obj.peakCounts = peakCounts;
            obj.spikeCounts = spikeCounts;
            
            % Resampled only.
            obj.signalRaw = signal;
            
            % Filtered.
            obj.fSmooth = fSmooth;
            
            obj.f0 = f0;
            obj.f1 = f1;
            obj.dff = dff;
            obj.spikes = spikes;
            obj.spikePeakIds = spikePeakIds;
            obj.area = area;
            obj.duration = duration;
            
            obj.nConditions = nConditions;
            obj.nCells = nCells;
            obj.cleanId = cleanId;
            obj.peakThreshold = peakThreshold;
            obj.peakMaskAll = peakMaskAll;
            obj.boxplotIds = boxplotIds;
            obj.boxplotLabels = boxplotLabels;
            obj.normalizedArea = normalizedArea;
            obj.firingRate = firingRate;
        end
        
        function plot(obj)
            obj.plotTrace();
            obj.plotTriggerAverage();
        end
        
        function fig = plotTrace(obj)
            % Plot traces and event windows.
            name = 'Traces and events';
            fig = figure('name', name);
            hold('all');
            % Add space between traces.
            pad = 5 * std(obj.dff(:));
            separation = pad * ones(1, obj.nCells);
            separation(1) = 0;
            separation = cumsum(separation);
            ylims = [separation(1) + min(obj.dff(:, 1)), separation(end) + max(obj.dff(:, end))];
            
            % Spikes.
            if ~isempty(obj.spikes)
                for s = 1:obj.nCells
                    xx = obj.time(find(obj.spikePeakIds(:, s)));
                    n = numel(xx);
                    xx = [xx, xx, NaN(n, 1)]';
                    yy = repmat([[s - 1, s] * pad - 0.5 * pad, NaN], n, 1)';
                    h = plot(xx(:), yy(:), 'LineStyle', '-', 'Color', 0.75 * ones(1, 3), 'Marker', 'none', 'DisplayName', 'Spikes', 'HandleVisibility', 'off');
                end
                h.HandleVisibility = 'on';
            end
            
            % Traces.
            for c = 1:obj.nConditions
                epochName = obj.configuration.conditionEpochs{2 * c - 1};
                epochTimes = obj.configuration.conditionEpochs{2 * c};
                epochTimes(epochTimes == -Inf) = min(obj.time);
                epochTimes(epochTimes == +Inf) = max(obj.time);                
                [faces, vertices] = patchEpochs(epochTimes, ylims(1), ylims(2));
                if isempty(obj.spikes)
                    label = sprintf('%s | peaks:%i', epochName, sum(obj.peakCounts(c, :)));
                else
                    label = sprintf('%s | peaks:%i | spikes:%i', epochName, sum(obj.peakCounts(c, :)), sum(obj.spikeCounts(c, :)));
                end
                patch('Faces', faces, 'Vertices', vertices, 'FaceColor', obj.cmap(c, :), 'EdgeColor', 'none', 'FaceAlpha', 0.50, 'DisplayName', label);
            end
            xx = repmat(obj.time, 1, obj.nCells);
            yy = bsxfun(@plus, separation, obj.dff);
            plot(xx, yy, 'LineWidth', 0.5, 'HandleVisibility', 'off');
            
            % Peaks.
            plot(xx(obj.peakMaskAll), yy(obj.peakMaskAll), 'k.', 'HandleVisibility', 'off');
            
            axis('tight');
            legend('show');
            set(gca, 'YTick', separation, 'YTickLabel', 1:obj.nCells);
            xlabel('Time (s)');
            ylabel('Cell number');
        end
        
        function figs = plotTriggerAverage(obj)
            figs(1) = figure('name', 'MSA: Peak-triggered average');
            [~, ~, patches] = triggerAverage(obj.dff, obj.peakIds, obj.peakLabels, obj.windowTemplate, obj.frequency, [], obj.configuration.conditionEpochs(1:2:end), obj.cmap);
            % Append width and height stats to peak-triggered average plot.
            for c = 1:obj.nConditions
                k = obj.peakLabels == c;
                data = obj.peakWidth(k(1:numel(c))); % !!
                widthMean = mean(data);
                n = numel(widthMean);
                widthSEM = std(data) / sqrt(n);
                patch = patches{c};
                append = sprintf(', width=%f±%f (s)', widthMean, widthSEM);
                patch.DisplayName = [patch.DisplayName, append];
            end
            ylabel('df/f');
            title('Peak-triggered average');
            
            if numel(obj.eventIds) > 0
                figs(1) = figure('name', 'MSA: Event-triggered average');
                triggerAverage(obj.dff, obj.eventIds, obj.eventLabels, obj.windowTemplate, obj.frequency, [], obj.configuration.conditionEpochs(1:2:end), obj.cmap);
                ylabel('df/f');
                title('Event-triggered average');
            end
        end

        function [av, sem] = eventTriggerAverage(obj, time)
            if numel(obj.eventIds) > 0
                [av, sem] = triggerAverage(obj.dff, obj.eventIds, obj.eventLabels, obj.windowTemplate, obj.frequency, time, obj.configuration.conditionEpochs(1:2:end));
            else
                error('No events have been defined.');
            end
        end

        function plotMUA(obj, binSize)
            binWidth = round(binSize * obj.frequency);
            % MUA at highest sampling rate.
            [nRows, nCols] = size(obj.peakMaskAll);
            id = zeros([nRows, nCols], 'uint64');
            id(obj.peakMaskAll) = find(obj.peakMaskAll);
            r = uint64((0:nCols - 1) * nRows);
            id = id - r;
            id = id(id > 0);
            obj.muaCounts = histcounts(id, 'BinWidth', binWidth);
            
            % Plot.
            figure('name', 'MSA: Multi-unit activity');
            
            % % Option 1 - line plot.
            % y = repmat(y, binWidth, 1);
            % y = y(1:nRows);
            % plot(obj.time, y);
            
            % Option 2 - bar plot.
            obj.muaTime = obj.time(1:binWidth:end);
            bar(obj.muaTime, obj.muaCounts, 1.0);
            
            xlabel('Time (s)');
            ylabel('MUA');
            axis('tight');
        end
        
        function export(obj, prefix)
            prefix = regexprep(prefix, '/$', '');
            [folder, basename] = fileparts(prefix);
            
            % Save peak counts per epoch and cell to file.
            output = fullfile(folder, sprintf('%s - peak count.csv', basename));
            fid = fopen(output, 'w');
            fprintf(fid, '# Peak count. Columns = cells. Rows = epochs.\n');
            fprintf(fid, [strjoin(repmat({'%d'}, 1, obj.nCells), ',') '\n'], obj.peakCounts');
            fclose(fid);
            
            % Save spike counts per epoch and cell to file.
            output = fullfile(folder, sprintf('%s - spike count.csv', basename));
            fid = fopen(output, 'w');
            fprintf(fid, '# Spike count. Columns = cells. Rows = epochs.\n');
            fprintf(fid, [strjoin(repmat({'%d'}, 1, obj.nCells), ',') '\n'], obj.spikeCounts');
            fclose(fid);
        end
    end
end

function output = parseNormalization(parameters, f, time)
    if iscell(parameters)
        fcn = parameters{1};
        if numel(parameters) == 1
            parameters{2} = [-Inf, Inf];
        end
        if isscalar(parameters{2})
            % Produce a vector from moving window.
            if numel(parameters) <= 2
                options = {'EndPoints', 'shrink'};
            else
                options = parameters(3:end);
            end
            frequency = 1 / median(diff(time));
            nSamples = numel(time);
            window = parameters{2};
            window = min(round(window * frequency), nSamples);
            output = fcn(f, window, options{:});
        else
            % Produce a value from all data (or epochs).
            epochs = parameters{2};
            ids = time2id(time, epochs);
            output = fcn(f(ids));
        end
    elseif isa(parameters, 'function_handle')
        % Produce a value from all data (or epochs).
        fcn = parameters;
        epochs = [-Inf, Inf];
        ids = time2id(time, epochs);
        output = fcn(f(ids));
    else
        output = parameters;
    end
end

function value = threshold(parameters, data)
    % {value1, @mad, @median}
    % {value1, @mad}
    % {value1, @mad, value2}
    % {value1}
    % value
    if iscell(parameters)
        n = numel(parameters);
        if n >= 1
            k = parameters{1};
        else
            k = 2.91;
        end
        if n >= 2
            f2 = parameters{2};
        else
            f2 = @mad;
        end
        if n >= 3
            f3 = parameters{3};
        else
            f3 = @median;
        end
    else
        k = parameters;
        f2 = @mad;
        f3 = @median;
    end
    if isa(f3, 'function_handle')
        value = k * f2(data) + f3(data);
    else
        value = k * f2(data) + f3;
    end
end

function output = warn(format, varargin)
    output = sprintf('[%s] %s', mfilename(), sprintf(format, varargin{:}));
end

function ylims = limits(x, percentile, grow)
    x = x(:);
    ylims = [prctile(x, 100 * (1 - percentile)), prctile(x, 100 * percentile)];
    delta = diff(ylims) * grow;
    ylims = [ylims(1) - delta, ylims(2) + delta];
end

function [odd, half] = forceOdd(fractional)
    % Number of samples in a triggered window (whole length, left to right).
    odd = fractional;
    % Force odd count.
    odd = round(odd) + (mod(round(odd), 2) == 0);
    half = (odd - 1) / 2;
end

function [averages, sems, patches] = triggerAverage(data, ids, labels, window, frequency, time, names, colors)
    % Filter out out-of-range traces.
    nSamples = numel(data);
    triggeredWindow = numel(window);
    halfWindow = (triggeredWindow - 1) / 2;
    k = ids > halfWindow & ids + halfWindow < nSamples;
    labels = labels(k);
    ids = ids(k);
    rawTime = window / frequency;
    nConditions = numel(names);
    if numel(time) > 0
        time = time(time >= rawTime(1) & time <= rawTime(end));
        averages = NaN(numel(time), nConditions);
    else
        averages = NaN(numel(rawTime), nConditions);
    end
    sems = averages;
    patches = cell(nConditions, 1);
    plotting = exist('colors', 'var');
    if numel(ids) > 0
        for c = 1:nConditions
            triggerIds = ids(labels == c);
            nTriggers = numel(triggerIds);
            if nTriggers > 0
                windowIds = triggerIds + window;
                triggeredData = data(windowIds);
                % Make sure matrix is nr x nc, particularly for 1 x nc.
                triggeredData = reshape(triggeredData, size(windowIds));
                nTriggers = size(triggeredData, 1);
                av = mean(triggeredData, 1);
                sem = std(triggeredData, [], 1) / sqrt(nTriggers);
                sem0 = sem(ceil(size(triggeredData, 2) / 2));
                if numel(time) > 0
                    av = interp1(rawTime, av, time);
                    sem = interp1(rawTime, sem, time);
                    t = time;
                else
                    t = rawTime;
                end
                averages(:, c) = av;
                sems(:, c) = sem;
                
                % Plot.
                if plotting
                    hold('all');
                    plot(t, av, 'Color', colors(c, :), 'HandleVisibility', 'off');
                    label = sprintf('%s n=%i, SEM=%f', names{c}, nTriggers, sem0);
                    vertices = [t; av + sem / 2];
                    vertices = cat(2, vertices, [fliplr(t); fliplr(av - sem / 2)])';
                    faces = 1:2 * triggeredWindow;
                    patches{c} = patch('Faces', faces, 'Vertices', vertices, 'FaceColor', colors(c, :), 'EdgeColor', 'none', 'FaceAlpha', 0.10, 'DisplayName', label);
                end
            end
        end
    elseif plotting
        text(0.5, 0.5, 'No triggers', 'HorizontalAlignment', 'center');
    end
    if plotting
        legend('show');
        xlabel('Time (s)');
        axis('tight');
    end
end