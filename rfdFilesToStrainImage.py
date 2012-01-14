#!/usr/bin/python
'''Take 3 filenames as input:
    pre RFD file
    post RFD file
    output strain file'''

from blockMatch import blockMatchClass
from sys import argv
import os

if len(argv) < 4:
    print "Error in usage.  Correct usage: "
    print os.path.basename(argv[0]) + " rfdFile1 rfdFile2 itkFile"

bl = blockMatchClass(argv[1],'rfd', argv[2] )
bl.CreateStrainImage(itkFileName=argv[3] )
