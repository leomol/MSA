% Load data; single unit activity recorded with confocal microscope and preprocessed with NeuroSeg.
inputDataFile = '../data/Confocal.csv';
data = loadData(inputDataFile);

% General configuration.
configuration = struct();
configuration.conditionEpochs = {'Baseline A', [0, 200], 'Drug A', [200, 250], 'Drug B', [400, Inf]};

% Peaks are 2 std from the mean of the first 100 seconds.
configuration.threshold = {2.00, @std, @mean};
configuration.thresholdEpochs = [0.00, 100.00];
configuration.triggeredWindow = 1.50;
configuration.peakDetectionMode = 'prominence';
configuration.lowpassFrequency = 5.0;
configuration.artifactEpochs = [400, 500];

% In the comments below:
%   f is the calcium response.
%   f0 is a value calculated from a subset of f.
%   window: length of the window enclosing neighbors of f.
normalizationOption = 5;
switch normalizationOption
    case 1
        % Do not normalize (default).
        configuration.f0 = 0;
        configuration.f1 = 1;
    case 2
        % z-score ==> (f - mean(f0)) / std(f0)
        configuration.f0 = @mean;
        configuration.f1 = @std;
    case 3
        % "modified z-score" ==> (f - median(f0)) / median(f0)
        configuration.f0 = @median;
        configuration.f1 = @mad;
    case 4
        % df/f ==> (f - f0) / f0
        configuration.f0 = @mean;
        configuration.f1 = @mean;
    case 5
        % z-score from a moving window of 60 seconds
        window = 60;
        configuration.f0 = {@movmean, window};
        configuration.f1 = {@movstd, window};
    case 6
        % df/f from a moving window of 60 seconds
        window = 60;
        configuration.f0 = {@movmean, window};
        configuration.f1 = {@movmean, window};
end

% Provide spikes or a function to deconvolve the flourescence from each cell.
if exist('deconvolveCa.m', 'file') == 2
    configuration.spikes = @deconvolve;
else
    configuration.spikes = [];
    warning('[%s] Could not find deconvolveCa.m', mfilename());
end

% Call MSA with given configuration.
msa = MSA(data, configuration);
% Print warnings, if any.
cellfun(@warning, msa.warnings);

% Plot trace window.
msa.plotTrace();
msa.plotTriggerAverage();

% Export data.
msa.export('exported');

function spikes = deconvolve(trace)
    [~, spikes, ~] = deconvolveCa(trace, 'foopsi', 'ar1', 'smin', -3, 'optimize_pars', true, 'optimize_b', true);
end