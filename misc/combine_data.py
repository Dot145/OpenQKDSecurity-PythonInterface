#!/usr/bin/python
import glob
import os
import numpy as np
from scipy.io import savemat

filenames = []
for name in glob.glob('*.dat'):
   filenames.append(name)

file = filenames[0]
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
      loss = 10**(loss/10) # and this converts from dB to actuall loss value
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
      print('Processed ' + file[:-4] + ".mat.")

savemat('alldata.mat', huge_dict)
print('Saved alldata.mat.')