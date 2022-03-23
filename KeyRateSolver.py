import numpy as np
import glob
import os
from scipy.io import savemat
import matlab.engine
import math

class KeyRateSolver:
    # constructor for KeyRateSolver object:
    # takes in a path, which should be a folder containing the .dat files that hold the data to be imported
    def __init__(self, path_to_data):
        if path_to_data[-1] == '/':
            self.path_to_data = path_to_data
        else:
            self.path_to_data = path_to_data + '/'
        # default time middle 347 (data row of t=0)
        self.time_middle = 347
        # default time range: 345-350
        self.time_range = np.arange(345, 351, 1, dtype='int')
        # default basis choice probabilities
        self.setReceiverBasisChoice('x', 0.666)

    def getKeyRate(self):
        if hasattr(self, 'eng'):
            self.createPreset()
            self.result = self.eng.getKeyRate46(self.path_to_data+'/alldata.mat')
            return self.result
        else:
            print('Error: please start the MATLAB engine by calling startEngine() on this object!')

    def startEngine(self):
        print('Starting MATLAB engine... (this may take a while)')
        self.eng = matlab.engine.start_matlab()
        self.eng.addpath(self.eng.genpath('code'))
        print('Done!')


    # function to let the user set the time range; we adjust by time_middle so that the user time is centered as
    # t=0 being the time of closest approach
    def setTimeRange(self, start, end, step):
        try:
            # clean up any decimals before passing into arange
            start = int(math.floor(start))
            step = int(round(step))
            # step of 0 doesn't work
            if step == 0:
                step = 1
            end = int(math.ceil(end))
            assert start < end
        except:
            print('Error: unable to set the time range. Make sure the starting time is less than the ending time, the step is positive, and that all arguments are integers.')
            print('\tUsage: setTimeRange(start, end, step)')
            return
        self.time_range = np.arange(start, end, step, dtype='int') + self.time_middle

    # function to let the user set the basis choice probabilities for Bob
    # the user choses a basis 'x', 'y', or 'z' as well as the probability p for that basis;
    # the other two basis choices will automatically be set to (1-p)/2
    # (any further asymmetry, though interesting, would not benefit key rate)
    def setReceiverBasisChoice(self, basis, p):
        if p > 0 and p < 1:
            if basis == 'x':
                self.pxB = p
                self.pzB = (1.-p)/2.
            elif basis == 'z':
                self.pzB = p
                self.pxB = (1.-p)/2.
            elif basis == 'y':
                self.pzB = (1.-p)/2.
                self.pxB = (1.-p)/2.
            else:
                print('Error: basis choice must be "x", "y", or "z"!')
        else:
            print('Error: primary basis choice probability must be between 0 and 1')
        

    # given a path to a folder containing a collection of .dat files, scans them and produces a matlab
    # .mat file that contains all of the info needed for our solver.
    def combineData(self, verbose = False):
        # get all .dat files in the directory
        filenames = []
        for name in glob.glob(self.path_to_data+'*.dat'):
            filenames.append(name)
        # check to make sure that there are .dat files, then parse them
        if len(filenames) > 0:
            huge_dict = {}
            for file in filenames:
                # count lines, but disregard the first two lines
                num_lines = sum(1 for line in open(file)) - 2
                with open(file) as f:
                    params = f.readline().split('\t')
                    paramDict = {}
                    for param in params:
                        key = param.split('=')[0]
                        try:
                            val = float(param.split('=')[1])
                        except:
                            val = param.split('=')[1]
                        paramDict[key] = val
                    # read the loss value in the file name
                    loss = file.split('_')[1] # this grabs "XdB" where X is the number we want
                    loss = float(loss[:-2]) # this grabs the number and makes it a float as well
                    loss = 1-10**(-loss/10) # and this converts from dB to actual loss value (1-effective transmittance)
                    # discard second line
                    f.readline()
                    times = np.zeros(num_lines)
                    elevations = np.zeros(num_lines)
                    distances = np.zeros(num_lines)
                    totChannels = np.zeros(num_lines)
                    detections = np.zeros((num_lines, 64))
                    i = 0
                    for line in f:
                    # line = f.readline()
                        tokens = line.split('\t')
                        [times[i], elevations[i], distances[i], totChannels[i]] = tokens[0:4]
                        detections[i,:] = tokens[4:]
                        i += 1
                    paramDict.update({'loss': loss, 'times': times, 'elevation': elevations, 'distance': distances, 'totChannel': totChannels, 'detections': detections})
                    #   savemat(file[:-4]+'.mat', paramDict)
                    bv = [paramDict['Bloch_x_coord'], paramDict['Bloch_y_coord'], paramDict['Bloch_z_coord']]
                    signal = 'unknown'
                    if bv == [1., 0., 0.]:
                        signal = 'D'
                    elif bv == [-1., 0., 0.]:
                        signal = 'A'
                    elif bv == [0., 1., 0.]:
                        signal = 'R'
                    elif bv == [0., -1., 0.]:
                        signal = 'L'
                    elif bv == [0., 0., 1.]:
                        signal = 'H'
                    elif bv == [0., 0., -1.]:
                        signal = 'V'
                    paramDict.update({'signal': signal})
                    key = 'mu_'+str(paramDict['mean_photon_no'])[-1]+'_'+signal
                    huge_dict[key] = paramDict
                    if verbose:
                        print('Processed ' + file[:-4] + ".mat.")

            savemat(self.path_to_data+'alldata.mat', huge_dict)
            print('Saved '+self.path_to_data+'alldata.mat.')
        else:
            print('Error: no .dat files found!')

    # print out appropriate .m file, including the times contained in the object
    def createPreset(self):
        # this is the first part of the preset file, up to declaring the time.
        rawPresetStr = '''
function [protocolDescription,channelModel,leakageEC,parameters,solverOptions] = SixStateDecoy46_asymptotic(decoys, etad, pzA, pzB, pxB, pd)
    [protocolDescription,channelModel,leakageEC]=setDescription();
    parameters=setParameters(decoys, etad, pzA, pzB, pxB, pd);
    solverOptions=setOptions();
end
function [protocolDescription,channelModel,leakageEC]=setDescription()
    
    protocolDescription=str2func('SixStateReducedLossyDescription');
    channelModel=str2func('SixStateDecoyChannel'); 
    leakageEC=str2func('generalEC');

end
function parameters=setParameters(decoys, etad, pzA, pzB, pxB, pd)

    parameters.names = ["etad", "pzB","pzA", "pxB","pd","decoys", "f", 'fullstat', 'time']; 
    parameters.scan.time = '''
        # now is where the time range is set
        rawPresetStr += np.array2string(self.time_range).replace('\n','')
        rawPresetStr += ''';
    parameters.fixed.pzA = pzA; 
    parameters.fixed.pzB = '''
        rawPresetStr += str(round(self.pzB, 4))
        rawPresetStr += ''';
    parameters.fixed.pxB = '''
        rawPresetStr += str(round(self.pxB, 4))
        # now print the rest of the file

        rawPresetStr += ''';
    parameters.fixed.pd = pd;
    parameters.fixed.f = 1.16;
    parameters.fixed.fullstat = 1;
    parameters.fixed.etad = etad;
    parameters.fixed.decoys = decoys;
end

function solverOptions=setOptions()
    solverOptions.globalSetting.cvxSolver = 'mosek';
    solverOptions.globalSetting.cvxPrecision = 'high';
    
    solverOptions.globalSetting.verboseLevel = 2; 
  
    solverOptions.optimizer.name = 'coordinateDescent'; 
    solverOptions.optimizer.linearResolution = 3; 
    solverOptions.optimizer.maxIterations = 2; 
    solverOptions.optimizer.linearSearchAlgorithm = 'iterative'; 
    solverOptions.optimizer.iterativeDepth = 2; 
    solverOptions.optimizer.maxSteps = 10; 
    solverOptions.optimizer.optimizerVerboseLevel = 0; 

    solverOptions.solver1.name = 'asymptotic_inequality';
    solverOptions.solver1.maxgap = 1e-6; 
    solverOptions.solver1.maxiter = 20;%10;
    solverOptions.solver1.initmethod = 1; %minimizes norm(rho0-rho) or -lambda_min(rho), use method 1 for finite size, 2 for asymptotic v1
    solverOptions.solver1.linearconstrainttolerance = 1e-8;
    solverOptions.solver1.linesearchprecision = 1e-20;
    solverOptions.solver1.linesearchminstep = 1e-3;
    solverOptions.solver1.maxgap_criteria = false; %true for testing gap, false for testing f1-f0
    solverOptions.solver1.removeLinearDependence = 'none'; %'qr' for QR decomposition, 'rref' for row echelon form, empty '' for not removing linear dependence
    solverOptions.solver2.name = 'asymptotic_inequality';
    solverOptions.solver2.epsilon = 0;
    solverOptions.solver2.epsilonprime = 1e-12;
end
'''
        # write this huge text to a file
        filename = 'code/SixStateDecoy46_asymptotic.m'
        f = open(filename, 'w')
        f.write(rawPresetStr)
        f.close()
        print('Wrote ' + filename + '\n with time values ' + str(self.time_range - self.time_middle) + '\n relative to closest approach.')
