% Load data; single unit activity recorded with nVista and preprocessed with min1pipe.
inputDataFile = '../data/miniscope.csv';
ttlFilename = '../data/miniscopeTTL.csv';
data = loadData(inputDataFile);

% Load TTL data.
[ttlRise, ttlFall] = loadInscopixTTL(ttlFilename);
% Trigger average at the center of each stimulation.
eventTimes = mean([ttlRise, ttlFall], 2);
conditionEpochs = [ttlRise, ttlFall]';
conditionEpochs = conditionEpochs(:)';

% General configuration.
configuration = struct();
configuration.conditionEpochs = {'TTL stimulation', conditionEpochs, 'TTL baseline', conditionEpochs - 4};
configuration.resamplingFrequency = 20.0;
configuration.eventTimes = eventTimes;

% Peaks are 2 std from the mean.
configuration.threshold = {2.00, @std, @mean};

% Call MSA with given configuration.
msa = MSA(data, configuration);
% Print warnings, if any.
cellfun(@warning, msa.warnings);

% Plot.
msa.plotTrace();