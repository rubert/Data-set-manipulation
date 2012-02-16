#!/usr/bin/python
'''Take 3 filenames as input:
    pre RFD file
    post RFD file
    output strain image'''

from blockMatch import blockMatchClass
from sys import argv
import os

if len(argv) < 3:
    print "Error in usage.  Correct usage: "
    print os.path.basename(argv[0]) + " rfdFile1 rfdFile2 pngFile <maxStrain>"

if len(argv) > 3:
    vMax = float( argv[3] )

else:
    vMax = .02

bl = blockMatchClass(argv[1],'rfd')
bl.SetRoiBoxSelect()
bl.CreateStrainImage(pngFileName=argv[2] + str(fNo).zfill(3) + '.png', preFrame = fNo , vMax = .02, skipFrame = 1)
