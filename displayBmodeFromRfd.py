#!/usr/bin/python

from rfData import rfClass
import os,sys


if len(sys.argv) < 2 or len(sys.argv) > 3:
    print "\nError, Usage is: "
    print os.path.basename(sys.argv[0]) + " rfFile <filetype> \n"
    sys.exit()

if len(sys.argv) == 2:
    rf = rfClass(sys.argv[1] , 'rfd')

if len(sys.argv) == 3:
    rf = rfClass(sys.argv[1], sys.argv[2])

rf.DisplayBmodeImage()
