import numpy as np
import glob, re
import os
from scipy.io import savemat
import matlab.engine
import math

class KeyRateSolver:
    # constructor for KeyRateSolver object:
    # takes in a path, which should be a folder containing the .dat files that hold the data to be imported
    def __init__(self, path_to_data):
        self.path_to_data = path_to_data
        # default time range: 345-350
        self.time_range = np.arange(345, 351, 1, dtype='int')
        # allowed protocol (in case we add more in the future)
        self.PROTOCOLS = {'pm46', 'pm44', 'pm46finite', 'pm44finite'}
        # corresponding preset files
        self.PRESET_FILES = {'pm46': 'pm46Decoy_asymptotic', 'pm44': 'pmBB84Decoy_asymptotic', 'pm46finite': 'pm46Decoy_finite', 'pm44finite': 'pmBB84Decoy_finite'}
        # default protocol: pm46, but can also do pm44
        self.protocol = 'pm46'
        self.setProtocol(self.protocol)
        # default basis choice probabilities
        self.setReceiverBasisChoice('x', 0.666)

    def getKeyRate(self):
        if hasattr(self, 'eng'):
            self.createPreset()
            self.result = self.eng.getKeyRate46(self.path_to_data)
            return self.result
        else:
            print('Error: please start the MATLAB engine by calling startEngine() on this object!')

    def startEngine(self):
        print('Starting MATLAB engine... (this may take a while)')
        self.eng = matlab.engine.start_matlab()
        self.eng.addpath(self.eng.genpath('code'))
        print('Done!')

    def setProtocol(self, protocol):
        if protocol in self.PROTOCOLS:
            self.protocol = protocol
            # modify getKeyRate file to select the correct preset
            lines = open('code/getKeyRate46.m', 'r').readlines()
            lineNum = 0
            for line in lines:
                if re.search('preset *= *', line):
                    lines[lineNum] = '    preset = "' + self.PRESET_FILES[protocol] + '";\n'
                lineNum += 1
            out = open('code/getKeyRate46.m', 'w')
            out.writelines(lines)
            out.close()
        else:
            print('Error: ' + str(protocol) + ' is not a valid protocol.')
            print('Available protocols: ' + str(self.PROTOCOLS))

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
        self.time_range = np.arange(start, end, step, dtype='int')

    # function to let the user set the basis choice probabilities for Bob
    # the user choses a basis 'x', 'y', or 'z' as well as the probability p for that basis;
    # the other two basis choices will automatically be set to (1-p)/2
    # (any further asymmetry, though interesting, would not benefit key rate)
    # has two modes, one for 46 and one for 44; if 'y' is chosen as the basis for 44, an error will occur.
    def setReceiverBasisChoice(self, basis, p):
        # verify probability
        if p > 0 and p < 1:
            # for 4-6 protocol:
            if self.protocol == 'pm46' or self.protocol == 'pm46finite':
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
            # for BB84
            elif self.protocol == 'pm44' or self.protocol == 'pm44finite':
                if basis == 'x':
                    self.pxB = p
                    self.pzB = 1-p
                elif basis == 'z':
                    self.pzB = p
                    self.pxB = 1-p
                else:
                    print('Error: basis choice must be "x" or "z"! (currently using ' + self.protocol + ' protocol)')
        else:
            print('Error: primary basis choice probability must be between 0 and 1 exclusive')
        

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

    def createPreset(self):
        if self.protocol == 'pm46':
            self.createPreset46()
        elif self.protocol == 'pm44':
            self.createPreset44()
        elif self.protocol == 'pm46finite':
            self.createPreset46finite()
        elif self.protocol == 'pm44finite':
            self.createPreset44finite()

    # print out appropriate .m file, including the times contained in the object
    def createPreset46(self):
        # file name to write to
        filename = self.PRESET_FILES[self.protocol] + '.m'
        # open the file for writing
        writefile = open('code/'+filename, 'w')
        # read the default preset file for 4-6
        with open('code/preset_46.m', 'r') as f:
            # read each line in the file
            lines = f.readlines()
            for line in lines:
                # replace the time, pzB, and pxB parameters with what the user has set
                if re.search('parameters.scan.time', line):
                    writefile.write('    parameters.scan.time = ' + np.array2string(self.time_range).replace('\n','') + ';\n')
                elif re.search('parameters.fixed.pzB', line):
                    writefile.write('    parameters.fixed.pzB = ' + str(round(self.pzB, 4)) + ';\n')
                elif re.search('parameters.fixed.pxB', line):
                    writefile.write('    parameters.fixed.pxB = ' + str(round(self.pxB, 4)) + ';\n')
                else:
                    writefile.write(line)
        # close file
        writefile.close()
        print('Wrote ' + filename + '\n with time values ' + str(self.time_range) + '.')

    def createPreset46finite(self):
        # file name to write to
        filename = self.PRESET_FILES[self.protocol] + '.m'
        # open the file for writing
        writefile = open('code/'+filename, 'w')
        # read the default preset file for 4-6
        with open('code/preset_46_finite.m', 'r') as f:
            # read each line in the file
            lines = f.readlines()
            for line in lines:
                # replace the time, pzB, and pxB parameters with what the user has set
                if re.search('parameters.scan.time', line):
                    writefile.write('    parameters.scan.time = ' + np.array2string(self.time_range).replace('\n','') + ';\n')
                elif re.search('parameters.fixed.pzB', line):
                    writefile.write('    parameters.fixed.pzB = ' + str(round(self.pzB, 4)) + ';\n')
                elif re.search('parameters.fixed.pxB', line):
                    writefile.write('    parameters.fixed.pxB = ' + str(round(self.pxB, 4)) + ';\n')
                else:
                    writefile.write(line)
        # close file
        writefile.close()
        print('Wrote ' + filename + '\n with time values ' + str(self.time_range) + '.')


    def createPreset44(self):
        # file name to write to
        filename = self.PRESET_FILES[self.protocol] + '.m'
        # open the file for writing
        writefile = open('code/'+filename, 'w')
        # read the default preset file for BB84
        with open('code/preset_44.m', 'r') as f:
            # read each line in the file
            lines = f.readlines()
            for line in lines:
                # replace time and pz with user's parameters
                if re.search('parameters.scan.time', line):
                    writefile.write('    parameters.scan.time = ' + np.array2string(self.time_range).replace('\n', '') + ';\n')
                elif re.search('parameters.fixed.pz', line):
                    writefile.write('    parameters.fixed.pz = ' + str(round(self.pzB, 4)) + ';\n')
                else:
                    writefile.write(line)
        # close file
        writefile.close()
        print('Wrote ' + filename + '\n with time values ' + str(self.time_range) + '.')

    def createPreset44finite(self):
        # file name to write to
        filename = self.PRESET_FILES[self.protocol] + '.m'
        # open the file for writing
        writefile = open('code/'+filename, 'w')
        # read the default preset file for BB84
        with open('code/preset_44_finite.m', 'r') as f:
            # read each line in the file
            lines = f.readlines()
            for line in lines:
                # replace time and pz with user's parameters
                if re.search('parameters.scan.time', line):
                    writefile.write('    parameters.scan.time = ' + np.array2string(self.time_range).replace('\n', '') + ';\n')
                elif re.search('parameters.fixed.pz', line):
                    writefile.write('    parameters.fixed.pz = ' + str(round(self.pzB, 4)) + ';\n')
                else:
                    writefile.write(line)
        # close file
        writefile.close()
        print('Wrote ' + filename + '\n with time values ' + str(self.time_range) + '.')

