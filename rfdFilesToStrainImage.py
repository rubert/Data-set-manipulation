#!/usr/bin/python
'''Take 3 filenames as input:
    pre RFD file
    post RFD file
    output strain file'''

from blockMatch import blockMatchClass
from sys import argv
import os

if len(argv) < 3:
    print "Error in usage.  Correct usage: "
    print os.path.basename(argv[0]) + " rfdFile1 itkFile"

bl = blockMatchClass(argv[1],'rfd')

for fNo in range(bl.nFrames - 1):
    bl.CreateStrainImage(pngFileName=argv[2] + str(fNo).zfill(3) + '.png', preFrame = fNo , vMax = .02, skipFrame = 1)
