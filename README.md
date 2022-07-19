# OpenQKDSecurity-PythonInterface
This repository provides files for interfacing experimental data with https://github.com/nlutkenhaus/openQKDsecurity. In particular, this interface allows for experimental expectation data from e.g. a satellite to be used as an input to calculate key rates for a 4-6 protocol or a BB84 with passive detection and decoy state analysis.

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
The experimental data must be in a particular format to be read correctly. The current format is a MATLAB data struct (`.mat`) which contains a struct for each signal and decoy intensity (for example, if sending horizontal, vertical, diagonal, and anti-diagonal signals with four possible decoy intensities, there will be a total of 16 sub-structs). Each of these structs contains parameters relevant to the experiment which are read by the interface to generate an accurate key rate. The exact list of parameters needed is as follows:
- `times`: the list of times at which data was collected (ex: [1,2,3,4,5] for 5 sets of expectation data)
- `detections`: the actual (16 or 64) x (# time steps) detection data, given as probabilities (16 for BB84, 64 for 4-6 protocol)
- `mean_photon_no`: the mean photon number of the signal (AKA signal intensity)
- `signal`: a character denoting the signal polarization ('H', 'V', 'D', or 'A')
- `N`: (in the case of a finite size protocol) the number of signals sent corresponding to this intensity and signal choice. 
The file `data/BB84_testdata.mat` is an example of a file in the appropriate format. Note that if the user desires to switch protocols using the MATLAB version of the interface, they need only edit the line `preset = ...` in `getKeyRate46.m` (instructions in the comments). In addition, a MATLAB user should directly enter the basis choice probabilities in the preset files `SixStateDecoy46_asymptotic.m` (for 4-6) or `pmBB84Decoy_asymptotic.m` (for BB84). 

## KeyRateSolver Class Methods
Descriptions of the methods of the KeyRateSolver class, as well as an example of use, follow.
| Function | Inputs | Outputs | Description |
| ----------- | ----------- | ----------- | ----------- |
| \_\_init\_\_ | `path_to_data` : string | KeyRateSolver object | Initializes a KeyRateSolver object and stores the path to observation data (as a .mat file) |
| startEngine | none | none | starts the MATLAB engine in python. Must be called before getKeyRate() |
| setProtocol | `protocol` : string | none | Allows the user to switch between the BB84 and 4-6 protocols. Currently only supports asymptotic versions of both protocols. Available strings for protocols are `'pm44'` for BB84 and `'pm46'` for 4-6. Must be called before getPreset() or getKeyRate(), as the preset file depends on the correct protocol. |
| setTimeRange | `start`, `end`, `step` : ints | none | Sets the time range of the key rate solver. Closest approach of the satellite is at t=0. Functions identical to numpy's `arange` function. |
| setReceiverBasisChoice | `basis` : string, `probability` : int | none | Sets the basis choice probabilities for the receiver Bob. The first argument is the primary (biased) basis choice, whose probability will be set to the second argument `probability`. The other two basis choices are implicitly set to divide the remaining probability in half in the 4-6 protocol; for BB84, the other basis choice will, obviously, be set to the remaining probability exactly. |
| combineData | `verbose` : boolean (optional) | none | Scans the folder set in the constructor for .dat files and compiles them into a file called `alldata.mat`, which is used by the key rate solver. Not necessary for key rate calculation. |
| createPreset | none | none | Sets up a preset file in the solver's software with the time range set by setTimeRange. This function must be called each time the time range is changed, but do note that it is called in getKeyRate() |
| getKeyRate | none | dict | creates a preset file, reads the data file set by the constructor, then computes the key rate for the provided time range. The return value is a MATLAB struct, which is converted into to a python dictionary in the KeyRateSolver class, and it can also be accessed from the `result` field of a KeyRateSolver object. This result dictionary contains all of the parameters of the experiment as well as the results. An example of how to extract key rate from the result is given below. |

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

## Potential Issues and Workarounds
### I'm having problems installing the MATLAB engine in python!
There isn't a lot I can do here, but there are several options to try. One option that worked for me was to create a conda environment with an appropriate version of python (3.9.7, for example) and perform the setup there.

### I'm getting an error that says "The global CVX solver selection cannot be changed while a model is being constructed."
There are two things that can cause this to happen. The most common reason is that the program was aborted mid-computation, and so the convex optimization module CVX was aborted mid-computation, which can cause issues. This error can be solved by calling `startEngine()` again. 
Another source of this issue is when the imported data leads to an infeasible semidefinite program in the key rate solver, which typically means that the detection data is nonphysical or there is a mismatch between the experimental parameters (such as basis choice probabilities) that generated the data and the parameters contained in the input data.

### I'm getting a key rate of 0 at every step and the program is telling me there is a step 2 solver exception.
This is a bug in the main openQKDsecurity software when dealing with the 4-6 protocol. I have submitted a pull request to the main software repository to fix it, but until it's fixed, you can replace openQKDsecurity/Solvers/Asymptotic_Solver/step2solverAsymptotic.m with [this version](https://github.com/Dot145/openQKDsecurity/blob/preRelease/Solvers/Asymptotic_Solver/step2SolverAsymptotic.m).

## Loose Ends

The `misc` folder contains a python script that scans a folder for `.dat` files from RefQ and combines them into a `.mat` file to be used as an input to `getKeyRate46.m`. This is the same functionality provided by the combineData function of the KeyRateSolver class, but in the form of a standalone script.

Email scott.johnstun@uwaterloo.ca for any questions about setting up the interface or data formatting.


