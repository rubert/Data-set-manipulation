#!/usr/bin/python

import sys
import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import hann
from scipy import optimize
from chirpz import chirpz

#Size of the ROI used to calculate GS
WINDOWX = 25. #mm
WINDOWX_PA = 100 #A-lines
WINDOWY = 12. #mm
GS_BW = 0.80 #MHz,(radius) corresponds to scatterer spacing of .385 microns
LOW_CUT = 1.1 #MHz, there seems to be a huge peak at low frequencies between 1.1 and 1.4 MHz
HIGH_CUT = 1.4

if len(sys.argv) < 2:
    print '\nError, usage is: \n generalizedSpectrumAvgValue.py DATA_DIR <PHASED ARRAY?>'
    sys.exit()

if len(sys.argv) == 3:
    phased = True

else:
    phased = False

DATA_DIR = sys.argv[1]

import os

if not os.path.isdir(DATA_DIR):
    print '\nError ' + DATA_DIR + ' is not a valid path.'
    sys.exit()

from rfData import rfClass

##Get an ROI from the first file in the directory so I know
##its size in pixels
print '\n ROI selection to determine ROI size in pixels'
tmpRf = rfClass(DATA_DIR + '/' + os.listdir(DATA_DIR)[0], 'rfd')

tmpRf.SetRoiFixedSize(WINDOWX, WINDOWY)    
region = tmpRf.data[tmpRf.roiY[0]:tmpRf.roiY[1], tmpRf.roiX[0]:tmpRf.roiX[1]]
roiPoints, roiLines = region.shape
psdPoints = roiPoints//2

##Set up CZT parameters for PSD calculation### 
lowCut = 0.0
highCut = 15.0

freqStep = (highCut - lowCut)/psdPoints
spectrumFreq = np.arange(0,psdPoints)*freqStep  + lowCut

fracUnitCircle = (highCut - lowCut)/(tmpRf.fs/10**6)

cztW = np.exp(1j* (-2*np.pi*fracUnitCircle)/psdPoints )
cztA = np.exp(1j* (2*np.pi*lowCut/(tmpRf.fs/10**6) ) )


#First loop through the directories and try to find an appropriate center frequency
#start by manually selecting ROIs at similar depths
#Compute the mean center frequency of all these ROIs
caStatistic = []
counter = 0
allFiles = os.listdir(DATA_DIR )
allFiles = [f for f in allFiles if 'rfd' in f]

for fName in allFiles:

    tmpRf = rfClass(DATA_DIR + '/' + fName, 'rfd')
    if not phased:
        tmpRf.SetRoiFixedSize(WINDOWX, WINDOWY)    
    else:
        tmpRf.SetRoiFixedSize(WINDOWX_PA, WINDOWY)


    ###Save the ROI coordinates to a numpy array so I don't have to 
    ###keep clicking the same pictures over and over
    np.save(DATA_DIR + '/roi' + str(counter), np.array( [tmpRf.roiY[0], tmpRf.roiY[1], tmpRf.roiX[0], tmpRf.roiX[1]] ) )
    
    tmpRf.SaveRoiImage(DATA_DIR + '/roi' + str(counter) + '.png' )
    tmpRf.ReadFrame()
        
    #PSD fit
    windowFuncPsd = hann(psdPoints+2)[1:-1].reshape(psdPoints,1)
    powerSpectrum = np.zeros( psdPoints )
        
    region = tmpRf.data[tmpRf.roiY[0]:tmpRf.roiY[1], tmpRf.roiX[0]:tmpRf.roiX[1]]
    maxDataWindow = region
    for r in range(3):
        dataWindow = maxDataWindow[psdPoints*r:psdPoints*r + psdPoints, :]*windowFuncPsd
        fourierData = np.zeros( dataWindow.shape )
        for l in range(maxDataWindow.shape[1]):
            fourierData[:,l] = abs(chirpz(dataWindow[:,l], cztA, cztW, psdPoints))
        powerSpectrum += fourierData.mean(axis = 1)

    powerSpectrum/= maxDataWindow.shape[1]*3 

    #now fit the spectrum to a gaussian
    errfunc = lambda param,x,y: y - np.exp(- (x - param[0])**2/(2*param[1]**2) )
    param0 = (5., 3.0)
    args = (spectrumFreq,powerSpectrum/powerSpectrum.max() )
    param, message = optimize.leastsq(errfunc, param0, args)
    mu = param[0]
    sigma = param[1]

    ##Set up CZT parameters for GS calculation### 
    lowCut = mu - GS_BW
    highCut = mu + GS_BW

    freqStep = (highCut - lowCut)/roiPoints
    spectrumFreqs = np.arange(0,roiPoints)*freqStep  + lowCut

    fracUnitCircle = (highCut - lowCut)/(tmpRf.fs/10**6)

    cztW = np.exp(1j* (-2*np.pi*fracUnitCircle)/roiPoints )
    cztA = np.exp(1j* (2*np.pi*lowCut/(tmpRf.fs/10**6) ) )

    ####
    ####Now that I have the ROI, compute the system-normalized
    ####generalized spectrum of each one
    ####
    import math
    gsSysNorm = np.zeros((roiPoints, roiPoints ) , np.complex128)

    windowFuncGS = hann(roiPoints+2)[1:-1].reshape(roiPoints,1)
    dataWindow = region.astype('double')
    dataWindow*=windowFuncGS
    tempPoints, tempLines = dataWindow.shape
    
    for l in range(tempLines):
        fourierData = chirpz(dataWindow[:,l], cztA, cztW, roiPoints)
        outerProd = np.outer(fourierData, fourierData.conjugate() )
        gsSysNorm += outerProd/abs(outerProd)

    gsSysNorm /= tempLines

    #Compute collapsed average
    CA = np.zeros( len(spectrumFreqs) )
    countCA = np.zeros( len(spectrumFreqs) )

    for f1 in range( len(spectrumFreqs) ):
        for f2 in range( len(spectrumFreqs) ):
                
            delta = abs(f2 - f1)
            CA[delta] += abs(gsSysNorm[f1,f2])
            countCA[delta] += 1

    CA /= countCA
    
    ##Work out cut-off frequencies for collapsed Average
    ##(3/4 window size) and .385 microns
    fLowPowerSpectrumRollOff = 1540./(2*.75*WINDOWY*10**-3)*10**-6
    startIndPowerSpectrumRollOff = int( round(fLowPowerSpectrumRollOff / freqStep) )
    maxStartInd = int( round(LOW_CUT/freqStep) )
    stopInd = int( round(HIGH_CUT/freqStep)  )

    meanValueCalculationIndexes = np.array( range(startIndPowerSpectrumRollOff,maxStartInd) + range(stopInd, tempPoints)
    )

    ###Add statistic to list
    caStatistic += [ CA[maxStartInd:stopInd].max()/CA[meanValueCalculationIndexes].mean() ]

    ##Plot collapsed average save plot
    plt.plot( np.arange(0, freqStep*len(spectrumFreqs), freqStep) , CA)
    plt.yticks( [0.2, 0.4, 0.6, 0.8, 1.0] )
    plt.xlabel('Frequency (MHz)' )
    plt.ylabel('CA Magnitude' )
    plt.savefig(DATA_DIR + '/collapsedAverage' + str(counter) + '.png' )
    plt.show()
    counter += 1
    print 'ROI ' + str(counter) +  ' of ' + str(len(allFiles))


#Write statistic to text file
caStatistic = np.array( caStatistic )
f = open( DATA_DIR + '/results.txt' , 'w')
f.write( 'The MAX/MEAN statistic mean value is: ' + str(caStatistic.mean() )  +'\n')
f.write( 'The MAX/MEAN standard deviation is: ' + str(caStatistic.std() ) + '\n' )
f.write( 'The number of ROIs showing a peak in the range of interest is: ' + str(len(caStatistic[caStatistic > 2] ) ) +
'\n')
f.write( 'The number of ROIs showing a peak in the range of interest as a percentage is: ' +
str(len(caStatistic[caStatistic > 2] ) /float(len(caStatistic) )) )
f.close()
