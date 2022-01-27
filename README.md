# OpenQKDSecurity-PythonInterface
This repository provides files for interfacing experimental data with https://github.com/nlutkenhaus/openQKDsecurity. In particular, this interface allows for experimental expectation data from e.g. a satellite to be used as an input to calculate key rates for a 4-6 protocol with passive detection and decoy.

## Setup
1. This interface requires openQKDsecurity (https://github.com/nlutkenhaus/openQKDsecurity) to be installed, so begin by following the steps there. 
2. (optional) If you wish to call openQKDsecurity from Python, you will need to install the MATLAB Engine for Python by following the steps at https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html. This is not necessary if you wish to use this interface directly in MATLAB.
3. Place the files in the `code` folder from this repository in the openQKDsecurity folder.
4. To run the code from MATLAB, open `Main46.m` and set the `filename` parameter to the location of experimental data. To run the code from Python, execute the following commands in a python prompt in the folder containing `getKeyRate46.m`:
```
import matlab.engine # import MATLAB python engine
eng = matlab.engine.start_matlab()
filename = 'path/to/data.mat'
result = eng.getKeyRate46(filename)
```
where `filename` is the path to the MATLAB data file which contains experimental data.

## Data Format
The experimental data must be in a particular format to be read correctly. Though this format is subject to change (will likely become a JSON file in the near future), the current format is a MATLAB data struct (`.mat`) which contains a struct for each signal and decoy intensity (for example, if sending horizontal, vertical, diagonal, and anti-diagonal signals with four possible decoy intensities, there will be a total of 16 sub-structs). Each of these structs contains parameters relevant to the experiment such as mean photon number, misalignment, dark count, beamsplitter transmissivity, etc., which are read by the interface to generate an accurate key rate. Each struct also contains a list of times (which must be 1 by $n$) and detections (which must be $n$ by 64) for each of the $n$ time steps in the experimental data. 

Email sjohnstu@uwaterloo.ca for any questions about setting up the interface or data formatting.
