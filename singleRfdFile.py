#!/usr/bin/python
'''This script is for creating strain images
when only a single RFD file was recorded while
an electrode was continuously moved

Usage:

python singleRfdFile.py RFD_NAME
'''


MAX_STRAIN = .02

import os, sys
START_DIR = os.getcwd()
RFD_NAME = sys.argv[1]
BASE_DIR = RFD_NAME[:-4] + 'dir'

from rfData import rfClass
rf = rfClass(RFD_NAME, 'rfd')

os.mkdir(BASE_DIR)

#B-mode images
os.mkdir(BASE_DIR + '/Bmode')
for frame in range(rf.nFrames):

    rf.SaveBModeImages( BASE_DIR + '/Bmode/f' + str(frame).zfill(3) , frame)    

#use M-player to create a movie
os.chdir(BASE_DIR + '/Bmode')
os.system('mencoder "mf://*.png" -mf fps=4 -o ../bMode.avi -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=1600')
os.chdir(START_DIR)

#Create ITK strain images
os.mkdir(BASE_DIR + '/strain')
from blockMatch import blockMatchClass
bl = blockMatchClass(RFD_NAME, 'rfd')
bl.WriteParameterFile(BASE_DIR +  '/blockMatchParams.txt' )

for frame in range(rf.nFrames - 1):
    bl.CreateStrainImage(preFrame = frame,pngFileName = BASE_DIR + '/strain/frame' + str(frame).zfill(3) + '.png' , itkFileName =
    BASE_DIR + '/strain/frame' + str(frame).zfill(3), vMax = .010)

#use M-player to create a movie
os.chdir(BASE_DIR +'/strain')
os.system('mencoder "mf://*.png" -mf fps=4 -o ../strain.avi -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=1600')
os.chdir(START_DIR)
