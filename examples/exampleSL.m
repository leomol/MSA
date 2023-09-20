% 2023-09-13. LM.
% 2023-09-13. Last modified.

%% Load data.
% Load data; single unit activity recorded with nVista and preprocessed with min1pipe.
% inputDataFile = 'C:\Users\szlsh\OneDrive\Desktop\Grin Lens\Day6 HFD_Feb5,2023_Inscopix_data_07-06-2023\D6 HFD Celltrace_1stDataTime.csv';
inputDataFile = 'data\Cell Trace Paradigm D0 30min.csv';
data = loadInscopix(inputDataFile);

%% Process.
binSize = 20;
events = 700;

% General configuration.
configuration = struct();
configuration.conditionEpochs = {'Before', [0, 600], 'After', [700, 760]};
configuration.resamplingFrequency = 10.0;
configuration.events = events;
configuration.triggeredWindow = 400;
% Peaks are 2 std from the mean.
configuration.threshold = {2.00, @std, @mean};

% Call MSA with given configuration.
msa = MSA(data, configuration);
% Print warnings, if any.
cellfun(@warning, msa.warnings);

%% Plot.
close all
msa.plotTrace();
msa.plotTriggerAverage();
msa.plotMUA(binSize);

% Plot MUA manually.
figure();
plot(msa.muaTime, msa.muaCounts);

%% Plot mean for each event and condition (a column per condition as defined in conditionEpochs).
step = 1;
time = -configuration.triggeredWindow(1) / 2:step:configuration.triggeredWindow(1) / 2;
av = msa.eventTriggerAverage(time);
figure();
hold('all');
nEpochs = numel(configuration.conditionEpochs) / 2;
for c = 1:nEpochs
    plot(time, av(:, c), 'DisplayName', sprintf(configuration.conditionEpochs{2 * c - 1}));
end
legend('show')