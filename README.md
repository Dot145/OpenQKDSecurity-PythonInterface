# OpenQKDSecurity-PythonInterface
This repository provides files for interfacing experimental data with https://github.com/nlutkenhaus/openQKDsecurity. In particular, this interface allows for experimental expectation data from e.g. a satellite to be used as an input to calculate key rates for a 4-6 protocol with passive detection and decoy.

## Setup
1. This interface requires openQKDsecurity (https://github.com/nlutkenhaus/openQKDsecurity) to be installed, so begin by following the steps there. 
2. (optional) If you wish to call openQKDsecurity from Python, you will need to install the MATLAB Engine for Python by following the steps at https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html. This is not necessary if you wish to use this interface directly in MATLAB.
3. Place the files in the `code` folder from this repository in the openQKDsecurity folder.
4. To run the code from MATLAB, open `Main46.m` and set the `filename` parameter to the location of experimental data, then simply run the program. 
5. To run the code from Python, create a KeyRateSolver object. The following is an example of how that might be done:
```python
from KeyRateSolver import KeyRateSolver
solver = KeyRateSolver('path/to/data/')
```
where `filename` is the path to the MATLAB data file which contains experimental data.

## Data Format
The experimental data must be in a particular format to be read correctly. Though this format is subject to change (will likely become a JSON file in the near future), the current format is a MATLAB data struct (`.mat`) which contains a struct for each signal and decoy intensity (for example, if sending horizontal, vertical, diagonal, and anti-diagonal signals with four possible decoy intensities, there will be a total of 16 sub-structs). Each of these structs contains parameters relevant to the experiment such as mean photon number, misalignment, dark count, beamsplitter transmissivity, etc., which are read by the interface to generate an accurate key rate. Each struct also contains a list of times (which must be 1 by n) and detections (which must be n by 64) for each of the n time steps in the experimental data. The file `data/alldata.mat` is an example of a file in the appropriate format; it can be re-created by calling the `combineData` function of the KeyRateSolver object on a folder containing RefQ observation data (.dat files).

## KeyRateSolver Class Methods
Descriptions of the methods of the KeyRateSolver class, as well as an example of use, follow.
| Function | Inputs | Outputs | Description |
| ----------- | ----------- | ----------- | ----------- |
| \_\_init\_\_ | path to data : string | KeyRateSolver object | Initializes a KeyRateSolver object and stores the path to observation data (the folder that contains .dat files or `alldata.mat` |
| startEngine | none | none | starts the MATLAB engine in python. Must be called before getKeyRate() |
| setTimeRange | start, end, step : ints | none | Sets the time range of the key rate solver. Closest approach of the satellite is at t=0. Functions identical to numpy's `arange` function. |
| combineData | verbose : boolean (optional) | none | Scans the folder set in the constructor for .dat files and compiles them into a file called `alldata.mat`, which is used by the key rate solver |
| createPreset | none | none | Sets up a preset file in the solver's software with the time range set by setTimeRange. This function must be called each time the time range is changed, but do note that it is called in getKeyRate() |
| getKeyRate | none | result : dict | creates a preset file, reads `alldata.mat` in the data path set by the constructor, then computes the key rate for the provided time range. The return value is a MATLAB struct, which is converted into to a python dict in the KeyRateSolver class. This contains all of the parameters of the experiment as well as the results. An example of how to extract key rate from the result is given below. |

Here is an example of how to set up the KeyRateSolver class in Python, where the RefQ .dat files are contained in the `data/` directory:
```python
# import the KeyRateSolver object from the file
from KeyRateSolver import KeyRateSolver
# initialize object and point to data directory
k = KeyRateSolver('data/')
# set time range to be [-5, -3, -1, 1, 3, 5]
k.setTimeRange(-5, 6, 2)  
# combine .dat files into alldata.mat (not necessary if alldata.mat already exists and is correct)
k.combineData()
# create the preset file (technically not necessary here, as the getKeyRate() call below calls this as well)
k.createPreset()
# start the MATLAB engine
k.startEngine()
# calculate key rate and store the results
results = k.getKeyRate()
```

After this, we can extract key rate information from `results` and plot it, for example with matplotlib:
```python
import matplotlib.pyplot as plt
import numpy as np
# extract the times that were scanned over
times = np.array(results['parameters']['scan']['time'][0])
# these times are given in raw row value (i.e. t=0 is 347), so shift back to be centered at t=0:
times = times - 347
# extract the lower and upper bounds of the key rate
# (this is pretty convoluted becasue of MATLAB's restrictions on the format of the data that can be returned to Python)
lowerBound = np.array([results['results'][time_point]['lowerBound'] for time_point in range(len(results['results']))])
upperBound = np.array([results['results'][time_point]['upperBound'] for time_point in range(len(results['results']))])
# plot the data
plt.plot(times, lowerBound, label='Lower Bound')
plt.plot(times, upperBound, label='Upper Bound')
plt.legend()
plt.xlabel('Time from closest approach (s)')
plt.ylabel('Asymptotic Key Rate')
plt.show()
```

## Loose Ends

At the moment, the object overrides basis choice probabilities from the imported data and defaults to px = 2/3, py = pz = 1/6. This matches the data that has been generated so far, but is subject to change in the future (eventually, this overriding functionality will be removed and the basis choice probabilities will be directly imported).

The `misc` folder contains a python script that scans a folder for `.dat` files from RefQ and combines them into a `.mat` file to be used as an input to `getKeyRate46.m`. This is the same functionality provided by the combineData function of the KeyRateSolver class, but in the form of a standalone script.

Email scott.johnstun@uwaterloo.ca for any questions about setting up the interface or data formatting.


