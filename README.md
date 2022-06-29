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

## More
TBD

## License
© 2019 [Leonardo Molina][Leonardo Molina]

This project is licensed under the [GNU GPLv3 License][License].

[Leonardo Molina]: https://github.com/leomol
[MATLAB]: https://www.mathworks.com/downloads/
[License]: LICENSE.md