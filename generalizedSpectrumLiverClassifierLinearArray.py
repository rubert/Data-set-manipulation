#!/usr/bin/python

import sys
import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import hann
from scipy import optimize
from chirpz import chirpz

#Size of the ROI used to calculate GS
WINDOWX = 10. #mm
WINDOWY = 8. #mm

if len(sys.argv) < 2:
    print '\nError, usage is: \n generalizedSpectrumLiverClassifierLinearArray.py DATA_DIR'
    sys.exit()

DATA_DIR = sys.argv[1]
LIVER_DIR = DATA_DIR + '/liver'
PHANTOM_DIR = DATA_DIR + '/phantom'
ROI_DIR = DATA_DIR + '/roi'

import os

if not os.path.isdir(DATA_DIR):
    print '\nError ' + DATA_DIR + ' is not a valid path.'
    sys.exit()

RESULT_DIR = DATA_DIR + '/results'
'''
if os.path.isdir(RESULT_DIR):
    os.rmdir(RESULT_DIR)
os.mkdir(RESULT_DIR)

if os.path.isdir(ROI_DIR):
    os.rmdir(ROI_DIR)
os.mkdir(ROI_DIR)
'''
from rfData import rfClass

##Get an ROI from the first file in the directory so I know
##its size in pixels
print '\n ROI selection to determine ROI size in pixels'
tmpRf = rfClass(LIVER_DIR + '/' + os.listdir(LIVER_DIR)[0], 'rfd')
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
'''
#First loop through the directories and try to find an appropriate center frequency
#start by manually selecting ROIs at similar depths
#Compute the mean center frequency of all these ROIs
runningTotalMu = 0.
runningTotalSigma = 0.
counter = 0

for fName in os.listdir(DATA_DIR + '/liver'):

    tmpRf = rfClass(LIVER_DIR + '/' + fName, 'rfd')
    tmpRf.SetRoiFixedSize(WINDOWX, WINDOWY)    
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

    runningTotalMu += mu
    runningTotalSigma += sigma
    counter += 1

avgMu = runningTotalMu/counter
avgSigma = runningTotalSigma/counter

print '\nInitial region selection complete: \n Mean center frequency of: ' + str(avgMu) + '\n Mean bandwidth of: ' + str(avgSigma)


##
##  Now pick ROIs in the liver
##
counter = 0
for fName in os.listdir(LIVER_DIR):
    while True:
        tmpRf = rfClass( LIVER_DIR + '/' + fName, 'rfd')
        tmpRf.SetRoiFixedSize(WINDOWX, WINDOWY)    

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

        print '\n Center frequency for this ROI\'s power spectrum is: ' + str(mu)
        print 'Average center frequency in trial selection was: ' + str(avgMu)
        
        print '\n Bandwidth for this ROI\'s power spectrum is: ' + str(sigma)
        print 'Average bandwidth in trial selection was: ' + str(avgSigma)

        answer = raw_input('Am I satisfied with this ROI? Type Y or y for yes')
        if answer == 'Y' or answer == 'y':
            counter += 1
            satisfiedWithAnswer = True
            np.save('results/liverRoi' + str(counter), region )
            tmpRf.SaveRoiImage(ROI_DIR + '/liverRoi' + str(counter) + '.png')
            
        answer = raw_input('Move on to next image?  Type Y or y for yes')
        if answer == 'Y' or answer == 'y':
            break        
        
##
##  Now pick ROIs in the TM phantom, can have more than 1 ROI per
##  image
counter = 0

for fName in os.listdir(PHANTOM_DIR):
    while True:
        tmpRf = rfClass(PHANTOM_DIR + '/' + fName, 'rfd')
        tmpRf.SetRoiFixedSize(WINDOWX, WINDOWY)    

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

        print '\n Center frequency for this ROI\'s power spectrum is: ' + str(mu)
        print 'Average center frequency in trial selection was: ' + str(avgMu)
        
        print '\n Bandwidth for this ROI\'s power spectrum is: ' + str(sigma)
        print 'Average bandwidth in trial selection was: ' + str(avgSigma)
        
        answer = raw_input('Am I satisfied with this ROI? Type Y or y for yes')
        if answer == 'Y' or answer == 'y':
            counter += 1
            satisfiedWithAnswer = True
            np.save('results/phantomRoi' + str(counter), region )
            tmpRf.SaveRoiImage(ROI_DIR + '/phantomRoi' + str(counter) + '.png')
            
        answer = raw_input('Move on to next image?  Type Y or y for yes')
        if answer == 'Y' or answer == 'y':
            break        
'''
avgMu = 4.73
avgSigma = .67
##Set up CZT parameters for GS calculation### 
lowCut = avgMu - 1.5
highCut = avgMu + 1.5

freqStep = (highCut - lowCut)/roiPoints
spectrumFreqs = np.arange(0,roiPoints)*freqStep  + lowCut

fracUnitCircle = (highCut - lowCut)/(tmpRf.fs/10**6)

cztW = np.exp(1j* (-2*np.pi*fracUnitCircle)/roiPoints )
cztA = np.exp(1j* (2*np.pi*lowCut/(tmpRf.fs/10**6) ) )

####
####Now that I have all the ROIs, compute the system-normalized
####generalized spectrum of each one
####
'''counter = 0
allFiles = os.listdir('results')
for fName in allFiles:

    tmpRoi = np.load('results' + '/' + fName)
   
    import math
    gsSysNorm = np.zeros((roiPoints, roiPoints ) , np.complex128)
    gsStdDev = np.zeros((roiPoints, roiPoints ) , np.complex128)

    windowFuncGS = hann(roiPoints+2)[1:-1].reshape(roiPoints,1)
    dataWindow = tmpRoi.astype('double')
    dataWindow*=windowFuncGS
    tempPoints, tempLines = dataWindow.shape
    
    for l in range(tempLines):
        fourierData = chirpz(dataWindow[:,l], cztA, cztW, roiPoints)
        outerProd = np.outer(fourierData, fourierData.conjugate() )
        gsSysNorm += outerProd/abs(outerProd)

    gsSysNorm /= tempLines

    #Compute the standard deviation of the GS
    #for l in range(tempLines):
    #    fourierData = chirpz(dataWindow[:,l], cztA, cztW, roiPoints)
    #    outerProd = np.outer(fourierData, fourierData.conjugate() )
    #    gsStdDev +=  abs(gsSysNorm - outerProd/abs(outerProd) )**2 

    #gsStdDev = np.sqrt(gsStdDev)/(tempLines - 1)
    #Scale the system-normalized GS by its standard deviation
    #gsSysNorm /= gsStdDev

    np.save('results' + '/' + fName[:-4] + 'GS', gsSysNorm)
    counter += 1
    print 'ROI ' + str(counter) +  ' of ' + str(len(allFiles))
'''
####
####  Now create a liver template  
####
meanLiver = np.zeros((roiPoints, roiPoints ) )
stdLiver = np.zeros((roiPoints, roiPoints ) )
liverFiles = [fName for fName in os.listdir('results') if 'liver' in fName]
liverFiles = [fName for fName in liverFiles if 'GS' in fName]
for fName in liverFiles:
    meanLiver += abs( np.load('results' + '/' + fName) )

meanLiver /= len(liverFiles)

for fName in liverFiles:
    stdLiver += (meanLiver - abs(  np.load('results' + '/' +fName) ) )**2

stdLiver /= (len(liverFiles) - 1)

####
####  Now create a phantom template  
####
meanPhantom = np.zeros((roiPoints, roiPoints ) )
stdPhantom = np.zeros((roiPoints, roiPoints ) )

phantomFiles = [fName for fName in os.listdir('results') if 'phantom' in fName]
phantomFiles = [fName for fName in phantomFiles if 'GS' in fName]

for fName in phantomFiles:
    meanPhantom += abs( np.load('results' + '/' + fName) )

meanPhantom /= len(phantomFiles)

for fName in phantomFiles:
    stdPhantom += (meanLiver - abs( np.load('results' + '/' + fName) ) )**2

stdPhantom = np.sqrt(stdPhantom)/(len(phantomFiles) - 1)


###Exclude spacings greater than the window length from the analysis
###
fishersDiscriminant = np.zeros( meanLiver.shape)
for f1 in range( roiPoints):
    for f2 in range( roiPoints):
        if f1 == f2:
            continue
        if 1540./(2.* abs(f1-f2)*freqStep*10**6)*10**3 < 8. :
            fishersDiscriminant[f1, f2] = (meanLiver[f1,f2] - meanPhantom[f1,f2])/(stdPhantom[f1,f2] + stdLiver[f1,f2] )

###
###  Now classify ROIs as either liver or phantom
###

plt.imshow( np.flipud(fishersDiscriminant)/fishersDiscriminant.max() , extent = [spectrumFreqs.min(),
spectrumFreqs.max(), spectrumFreqs.min(), spectrumFreqs.max()])
plt.colorbar()
plt.savefig('results/linearDiscriminantImage.png')
plt.close()

phantomScores = np.zeros( len(phantomFiles) )
liverScores = np.zeros( len(liverFiles) )

for ind,fName in enumerate(phantomFiles):
    phantomScores[ind] = (fishersDiscriminant*abs(np.load('results' + '/' + fName))).sum()
    
for ind, fName in enumerate(liverFiles):
    liverScores[ind] = (fishersDiscriminant*abs(np.load('results' + '/' + fName))).sum()

phantomScores /= liverScores.max()
liverScores /= liverScores.max()

plt.plot(np.arange(1,len(phantomScores) + 1), phantomScores, 'ro', markersize = 12, label = 'Liver')
plt.plot(np.arange(1,len(liverScores) + 1), liverScores, 'bo', markersize = 12 ,label = 'Ablation')
plt.xlabel('Index of ROI')
plt.ylabel('Magnitude of template response')
plt.legend()
plt.savefig('results/score.png')

q = 0


