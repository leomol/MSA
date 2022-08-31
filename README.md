
# Single unit analisis
Script to process single unit data and quantify activity (e.g. spontaneous peaks and correlations) on arbitrary epoch definitions.
![MSA example](etc/example.png)

## Installation
- Download and install [MATLAB][MATLAB].
- Download and extract source files to `Documents/MATLAB`
- For analyses involving deconvolution of CA²⁺ data, download and implement such library in the code. For example, download [OASIS library](https://github.com/zhoupc/OASIS_matlab), run `oasis_setup.m`, and use one of the examples under [examples](examples).

## Usage
- Run `setup.m` to add FPA dependencies to the MATLAB search path.
- Create a copy of the example script.
- Change the configuration variables and the data loader according to your data acquisition system.
- Define epochs.
- Change settings and thresholds to compute peaks and clusters.
- Click `Run`.

Epochs consists of an array of timestamps indicating start and stop times (in seconds) of a given condition `{'epoch1', [start1, end1, start2, end2, ...], 'epoch2', ...}`. For example:
```MATLAB
  epochs = {'A', [100, 200], 'B', [200, 500], 'C', [500, 600, 600, 700]};
```
means there are three epochs called A, B, and C:
```
A: 100 to 200
B: 200 to 500
C: 500 to 600 then 600 to 700
```

## Documentation

`msa = MSA(data, configuration);`

Normalize, filter, and detect peaks of spontaneous activity in user defined epochs.

Each column in `data` corresponds to activity from individual cells.

`configuration` is a struct with the following fields (defaults are used for missing fields):
 - `conditionEpochs` - Epochs for different conditions: {'epoch1', [start1, end1, start2, end2, ...], 'epoch2', ...}
 - `thresholdEpochs` - Epochs to include for peak threshold calculation.
 - `events` - Event-triggered data; times at which a type of event occurs.
 - `resamplingFrequency` - Resampling frequency (Hz).
 - `lowpassFrequency` - Lowest frequency permitted in normalized signal.
 - `peakSeparation` - Minimum time separation between two peaks.

### Processing steps
- Resample signal and reference to a given frequency.
- Baseline correction modeled as an exponential decay of the low-pass filtered data (optionally using airPLS).
- Correct for motion artifacts by subtracting reference to signal, after a polynomial fit (optional).
- Remove fast oscillations with a low-pass filter.
- Normalize data as df/f or z-score according to settings.
- Find peaks of spontaneous activity in low-pass filtered data.

Normalization is calculated as (f - f0) / f1 where f0 and f1 can be data provided by the
user or calculated using given functions:

See examples and source code for detailed analysis steps and default parameters.
Units for time and frequency are seconds and hertz respectively.

## License
© 2019 [Leonardo Molina][Leonardo Molina]

This project is licensed under the [GNU GPLv3 License][License].

[Leonardo Molina]: https://github.com/leomol
[MATLAB]: https://www.mathworks.com/downloads/
[License]: LICENSE.md
