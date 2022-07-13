% Load data; single unit activity recorded with confocal microscope and preprocessed with NeuroSeg.
inputDataFile = '../data/Confocal.csv';
data = loadData(inputDataFile);

% General configuration.
configuration = struct();
configuration.conditionEpochs = {'Baseline A', [0, 200], 'Drug A', [200, 250], 'Drug B', [400, Inf]};

% Peaks are 2 std from the mean of the first 100 seconds.
configuration.threshold = {2.00, @std, @mean};
configuration.triggeredWindow = 1.50;
configuration.thresholdEpochs = [0.00, 100.00];
configuration.peakDetectionMode = 'prominence';

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