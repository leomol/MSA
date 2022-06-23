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
            defaults.conditionEpochs = {'Data', [-Inf, Inf]};
            defaults.artifactEpochs = [];
            defaults.eventTimes = [];
            defaults.resamplingFrequency = NaN;
            defaults.lowpassFrequency = Inf;
            defaults.peakSeparation = 2.0;
            defaults.f0 = 0.0;
            defaults.f1 = 1.0;
            defaults.thresholdEpochs = NaN;
            defaults.threshold = {2.91, @mad, @median};
            defaults.triggeredWindow = 10.0;
            defaults.spikes = [];
            
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
                if any(+dff(:, u) >= peakThreshold(u))
                    % !!
                    %[~, k, peakWidth] = findpeaks(+dff(:, u), time, 'MinPeakHeight', peakThreshold(u), 'MinPeakDistance', configuration.peakSeparation, 'WidthReference', 'halfheight');
                    %k = intersect(find(ismember(time, k)), thresholdId);
                    [~, k, peakWidthCell] = findpeaks(+dff(:, u), 'MinPeakHeight', peakThreshold(u), 'MinPeakDistance', configuration.peakSeparation * frequency, 'WidthReference', 'halfheight');
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
            x = arrayfun(@(t) find(time >= t, 1, 'first'), configuration.eventTimes, 'UniformOutput', false);
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
                [faces, vertices] = patchEpochs(epochTimes, ylims(1), ylims(2));
                if isempty(obj.spikes)
                    label = sprintf('%s | peaks:%i', epochName, sum(obj.peakCounts(c, :)));
                else
                    label = sprintf('%s | peaks:%i | firing rate:%.2f (Hz)', epochName, sum(obj.peakCounts(c, :)), mean(obj.firingRate(c, :)));
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
            patches = plotTriggerAverage(obj.dff, obj.peakIds, obj.peakLabels, obj.windowTemplate, obj.frequency, obj.configuration.conditionEpochs(1:2:end), obj.cmap);
            % Append width and height stats to peak-triggered average plot.
            for c = 1:obj.nConditions
                k = obj.peakLabels == c;
                data = obj.peakWidth(k);
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
                plotTriggerAverage(obj.dff, obj.eventIds, obj.eventLabels, obj.windowTemplate, obj.frequency, obj.configuration.conditionEpochs(1:2:end), obj.cmap);
                ylabel('df/f');
                title('Event-triggered average');
            end
            
            figs(2) = figure('name', 'MSA: start-triggered average');
            plotTriggerAverage(obj.dff, obj.epochIds(1:2:end)', obj.epochLabels, obj.windowTemplate, obj.frequency, obj.configuration.conditionEpochs(1:2:end), obj.cmap);
            ylabel('df/f');
            title('Start-triggered average');
            
            figs(3) = figure('name', 'MSA: stop-triggered average');
            plotTriggerAverage(obj.dff, obj.epochIds(2:2:end)', obj.epochLabels, obj.windowTemplate, obj.frequency, obj.configuration.conditionEpochs(1:2:end), obj.cmap);
            ylabel('df/f');
            title('Stop-triggered average');
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

function patches = plotTriggerAverage(data, ids, labels, window, frequency, names, colors)
    % Filter out out-of-range traces.
    nSamples = numel(data);
    triggeredWindow = numel(window);
    halfWindow = (triggeredWindow - 1) / 2;
    k = ids > halfWindow & ids + halfWindow < nSamples;
    labels = labels(k);
    ids = ids(k);
    time = window / frequency;
    nConditions = numel(names);
    patches = cell(nConditions, 1);
    if numel(ids) > 0
        hold('all');
        for c = 1:nConditions
            triggerIds = ids(labels == c);
            nTriggers = numel(triggerIds);
            if nTriggers > 0
                windowIds = triggerIds + window;
                triggeredData = data(windowIds);
                % Make sure matrix is nr x nc, particularly for 1 x nc.
                triggeredData = reshape(triggeredData, size(windowIds));
                av = mean(triggeredData, 1);
                nTriggers = size(triggeredData, 1);
                % Plot.
                plot(time, av, 'Color', colors(c, :), 'HandleVisibility', 'off');
                sem = std(triggeredData, [], 1) / sqrt(nTriggers);
                sem0 = sem(ceil(size(triggeredData, 2) / 2));
                label = sprintf('%s n=%i, SEM=%f', names{c}, nTriggers, sem0);
                vertices = [time; av + sem / 2];
                vertices = cat(2, vertices, [fliplr(time); fliplr(av - sem / 2)])';
                faces = 1:2 * triggeredWindow;
                patches{c} = patch('Faces', faces, 'Vertices', vertices, 'FaceColor', colors(c, :), 'EdgeColor', 'none', 'FaceAlpha', 0.10, 'DisplayName', label);
            end
        end
    else
        text(0.5, 0.5, 'No triggers', 'HorizontalAlignment', 'center');
    end
    legend('show');
    xlabel('Time (s)');
    axis('tight');
end