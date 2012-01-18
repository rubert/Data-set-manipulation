#!/usr/bin/python
'''  This script takes an RFD file as input.  It will create a processing directory, then create
B-mode images.  Two movies will be created, one at the original frame rate, one
at a frame rate of 5 Hz.'''

import sys, os

#Check input argument
if len(sys.argv) != 2:
    print "\nError, usage is: \n processByDQM.py RFD_NAME \n"
    sys.exit()

if not os.path.isfile(sys.argv[1]):
    print "\nError, usage is: \n processByDQM.py RFD_NAME "
    print "RFD_NAME passed is not a regular file.(Potential RFD file) \n"
    sys.exit()

#Establish paths
BASE_PATH, RFD_NAME = os.path.split(os.path.abspath(sys.argv[1]))
print "\nProcessing file named: " + RFD_NAME + "\non path: " + BASE_PATH + "\n"
RESULT_DIR = BASE_PATH +'/' + RFD_NAME[:-4] + 'Bmode'

if os.path.isdir(RESULT_DIR):
    print "\nError, data set already processed."
    print "Directory: " + RESULT_DIR + " Exists. \n"
    sys.exit()

os.makedirs(RESULT_DIR)
os.makedirs(RESULT_DIR + '/images')

from rfData import rfClass
rf = rfClass(BASE_PATH + '/' + RFD_NAME, 'rfd')


for frames in range(rf.nFrames):
    rf.SaveBModeImage(RESULT_DIR + '/images/frame_' + str(frames).zfill(3), image = frames)

#Create movies
import readrfd
header = readrfd.headerInfo(BASE_PATH + '/' + RFD_NAME)
fps = int( header[3] )
import os
os.system('''mencoder "mf://''' + RESULT_DIR + '''/images/*.png" -mf type=png:fps=5 -ovc lavc -o '''+ RESULT_DIR + '/bMode_fps5.avi')
os.system('''mencoder "mf://''' + RESULT_DIR + '''/images/*.png" -mf type=png:fps=''' + str(fps)+''' -ovc lavc -o '''+ RESULT_DIR + '/bMode_fpsActual'
+ str(fps) +'.avi' )
