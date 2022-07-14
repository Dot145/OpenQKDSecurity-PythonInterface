import re

writefile = open('testout.m', 'w')
with open('SixStateDecoy46_asymptotic.m') as f:
    lines = f.readlines()
    for line in lines:
        if re.search('parameters.scan.time', line):
            writefile.write('\tparameters.scan.time = ' + np.array2string(self.time_range).replace('\n','') + ';\n')
        elif re.search('parameters.fixed.pzB', line):
            writefile.write('\tparameters.fixed.pzB = ' + str(round(self.pzB, 4)) + ';\n')
        elif re.search('parameters.fixed.pxB', line):
            writefile.write('\tparameters.fixed.pxB = ' + str(round(self.pxB, 4)) + ';\n')
        else:
            writefile.write(line)
