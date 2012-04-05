#!/usr/bin/python
'''Take 3 filenames as input:
    pre RFD file
    post RFD file
    output strain image
    It's a script for my attenuation estimation experiment to
    spit out a strain image given 1-frame pre and post compression
    RFD files.'''

from blockMatch import blockMatchClass
from sys import argv
import os
from matplotlib import pyplot as plt

if len(argv) < 3:
    print "Error in usage.  Correct usage: "
    print os.path.basename(argv[0]) + " rfdFile1 rfdFile2 pngFile <maxStrain>"

if len(argv) > 4:
    vMax = float( argv[4] )

else:
    vMax = .02

if '.png' not in argv[3]:
    picFName = argv[3] + '.png'
else:
    picFName = argv[3]

bl = blockMatchClass(argv[1],'rfd',postFile = argv[2])#, selectRoi = True)
bl.CreateStrainImage(pngFileName=picFName, vMax = vMax)

bl.DisplayScanConvertedStrainImage(vMax)

answer = None
while answer is not 'Y' and answer is not 'N':
    answer = raw_input('Is the plotting scale satisfactory? (Y or N)?\n')

if answer is 'Y':
    satisfiedWithScale = True
else:
    satisfiedWithScale = False

while not satisfiedWithScale:
    newMax = float( raw_input('Enter a new maximum strain value to display\n') )
    print '\n The new maximum is: ' + str(newMax)
    bl.RescaleRgbStrainImage(newMax)
    bl.DisplayScanConvertedStrainImage(newMax)
   
    answer = None
    while answer is not 'Y' and answer is not 'N':
        answer = raw_input('Is the plotting scale satisfactory? (Y or N)? \n')

    if answer is 'Y':
        satisfiedWithScale = True
    else:
        satisfiedWithScale = False
