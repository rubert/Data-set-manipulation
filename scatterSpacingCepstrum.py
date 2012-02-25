#!/usr/bin/python

import sys
import numpy as np
from matplotlib import pyplot as plt
from mtspec import mtspec
from chirpz import chirpz

#Size of the ROI used to calculate GS
WINDOWX = 10. #mm
WINDOWX_PA = 80 #A-lines
WINDOWY = 12. #mm
GS_BW = 0.80 #MHz,(radius) corresponds to scatterer spacing of .385 microns
LOW_CUT2 = .9 #MHz, there seems to be a huge peak at low frequencies between 


if len(sys.argv) < 3:
    print '\nError, usage is: \n generalizedSpectrumAvgValue.py DATA_DIR REFPHANFILE <PHASED ARRAY?>'
    sys.exit()

if len(sys.argv) == 4:
    phased = True

else:
    phased = False

DATA_DIR = sys.argv[1]
REFPHANFILE = sys.argv[2]

import os

if not os.path.isdir(DATA_DIR):
    print '\nError ' + DATA_DIR + ' is not a valid path.'
    sys.exit()

from rfData import rfClass

##Get an ROI from the first file in the directory so I know
##its size in pixels
tmpRf = rfClass(DATA_DIR + '/' + os.listdir(DATA_DIR)[0], 'rfd')

counter = 0
allFiles = os.listdir(DATA_DIR )
allFiles = [f for f in allFiles if 'rfd' in f]

refPhanRf = rfClass(REFPHANFILE, 'rfd')

for fName in allFiles:

    tmpRf = rfClass(DATA_DIR + '/' + fName, 'rfd')
    if not phased:
        tmpRf.SetRoiFixedSize(WINDOWX, WINDOWY)    
    else:
        tmpRf.SetRoiFixedSize(WINDOWX_PA, WINDOWY)
        refPhanRf.SetRoiFixedSize(WINDOWX_PA, WINDOWY)

    tmpRf.SaveRoiImage(DATA_DIR + '/roi' + str(counter) + '.png' )
    tmpRf.ReadFrame()
    
    refPhanRf.SaveRoiImage(DATA_DIR + '/Refroi' + str(counter) + '.png' )
    refPhanRf.ReadFrame()
        

    region = tmpRf.data[tmpRf.roiY[0]:tmpRf.roiY[1], tmpRf.roiX[0]:tmpRf.roiX[1]]
    refRegion = refPhanRf.data[refPhanRf.roiY[0]:refPhanRf.roiY[1], refPhanRf.roiX[0]:refPhanRf.roiX[1]]

    ####
    ####Now that I have the ROI, compute the power spectrum of the signal for
    ####each A-line using MTM and average
    ####
    import math
    dataWindow = region.astype('double')
    tempPoints, tempLines = dataWindow.shape
    halfPoints = int(tempPoints/2)
    if halfPoints%2:
        halfPoints += 1

    spectrum = np.zeros( halfPoints )
    
    for l in range(tempLines):
        tmpSpec, freq = mtspec(dataWindow[:,l],1./tmpRf.fs , 3.)  #data, delta, time-bandwidth
        spectrum += tmpSpec

    spectrum /= tempLines
    
    spectrumRef = np.zeros( halfPoints )
    dataWindowRef = region.astype('double')
    for l in range(tempLines):
        tmpSpec, freq = mtspec(dataWindowRef[:,l],1./tmpRf.fs , 3.)  #data, delta, time-bandwidth
        spectrumRef += tmpSpec

    spectrumRef /= tempLines

    ###
    maxF = freq[-1]
    deltaF = freq[1] - freq[0]
    deltaT = 1./(halfPoints*deltaF)

    ###Work out times I am interested in (cepstrum frequencies)
    highT = 2*1.00*WINDOWY*10**-3/1540.  ## windowY is in mm, use 3/4 window size for high cutoff
    lowT = 0#2*.4*10**-3/1540.  ##.4 mm low cutoff

    fracUnitCircle =(highT - lowT)/(deltaT*halfPoints)
    cztA = np.exp(1j* (2*np.pi*lowT/deltaT ) )
    cztW = np.exp(1j* (-2*np.pi*fracUnitCircle) )
   
    chirpTimes = np.linspace(lowT, highT, halfPoints)

    ##Compute Cepstrum
    #cepstrum = abs(chirpz( np.log( spectrum ), cztA, cztW , halfPoints))
    cepstrum, melFreq = mtspec(np.log(spectrum),deltaF, 3.)  #data, delta, time-bandwidth
    cepstrumRef, melFreq = mtspec(np.log(spectrumRef),deltaF, 3.)  #data, delta, time-bandwidth
    
    ##Plot collapsed average save plot
    plt.plot(melFreq*10.**6, cepstrum )
    plt.xlabel('Quefrency ($\mu$s)' )
    plt.savefig(DATA_DIR + '/cepstrum' + str(counter) + '.png' )
    plt.show()
    
    ##Plot collapsed average save plot
    plt.plot(melFreq*10.**6, cepstrum - cepstrumRef )
    plt.xlabel('Quefrency ($\mu$s)' )
    plt.savefig(DATA_DIR + '/cepstrum' + str(counter) + '.png' )
    plt.show()
    counter += 1
    print 'ROI ' + str(counter) +  ' of ' + str(len(allFiles))


    ###Subtract off the cepstrum from a similar ROI in a reference phantom
    


f = open( DATA_DIR + '/results.txt' , 'w')
f.write('Average CA of: ' + str(avgCA.mean()) +'\n' )
f.write( 'Standard dev. of: ' + str(avgCA.std())  )
f.write('\n\n')

f.write( 'In limited range, Average CA of: ' + str(avgCA2.mean())  +'\n' )
f.write( 'In limited range, Standard dev. of: ' + str(avgCA2.std())   )
f.write('\n\n')

f.write( 'Between ' + str(LOW_CUT2) + ' and 1.5 MHz the mean maximum value of the CA is: ' + str( np.array(caMax).mean() )  +'\n')
f.write( 'Between ' + str(LOW_CUT2) + ' and 1.5 MHz the std dev maximum value of the CA is: ' + str( np.array(caMax).std() )  )
f.write('\n\n')
f.close()
