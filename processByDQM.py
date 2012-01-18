#!/usr/bin/python
'''  This script takes an RFD file as input.  It will create a processing directory, then create
strain images, and select the best one based on the DQM.  It assumes that the user 
already knows the ROI they want from looking at the B-mode images.'''

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
RESULT_DIR = BASE_PATH +'/' + RFD_NAME[:-4] + 'Results'
if os.path.isdir(RESULT_DIR):
    print "\nError, data set already processed."
    print "Directory: " + RESULT_DIR + " Exists. \n"
    sys.exit()

os.makedirs(RESULT_DIR)
os.makedirs(RESULT_DIR + '/blockMatch')

#Perform block-matching
from blockMatch import blockMatchClass
bl = blockMatchClass(BASE_PATH +'/'+ RFD_NAME, 'rfd', selectRoi = True)

for frameNo in range(bl.nFrames):
    for skipNo in range(5):

        bl.CreateStrainImage(preFrame = frameNo, skipFrame = skipNo, itkFileName = RESULT_DIR + '/blockMatch/frame_' + str(frameNo) +'_'
        + str(frameNo + 1 + skipNo) )

#After block-matching I compute the frame triplets with the optimum DQM
import numpy
from numpy import zeros
import SimpleITK as sitk
import cv
from scipy.interpolate import interp1d, RectBivariateSpline
from interp2NoSplines import interp2
DQM = zeros( bl.nFrames )
framePairs = zeros( (bl.nFrames, 2) )

for frameNo in range(5, bl.nFrames -5):

    tmpDQM = zeros( (5,5) )
    
    rhoRfJ = zeros( (5,1) )
    rhoRfK = zeros( (1,5) )
    rhoRf = zeros( (5,5) )
    
    reader = sitk.ImageFileReader()
    
    for skip in range(5):
        #load up displacement between frames
        reader.SetFileName(RESULT_DIR + '/blockMatch/frame_' + str(frameNo) + '_' + str(frameNo + skip + 1) + 'dispY.mhd')
        dpYItk = reader.Execute()
        dpY = sitk.GetArrayFromImage(dpYItk)

        reader.SetFileName(RESULT_DIR + '/blockMatch/frame_' + str(frameNo) + '_' + str(frameNo + skip + 1) + 'dispX.mhd')
        dpXItk = reader.Execute()
        dpX = sitk.GetArrayFromImage(dpXItk)
       
        dpYUp = zeros( (bl.stopY - bl.stepY - bl.startY, bl.stopX - bl.startX) )
        dpXUp = zeros( (bl.stopY - bl.stepY - bl.startY, bl.stopX - bl.startX) )

        #Upsample displacement to same spacing as RF.
        paramRfYIndexes = numpy.arange(bl.startY, bl.stopY, bl.stepY )
        paramRfYIndexesNew = numpy.arange(bl.startY, bl.stopY - bl.stepY )

        #dpY
        for x in range(bl.stopX - bl.startX):
            interp = interp1d(paramRfYIndexes, dpY[:,x] , kind = 'cubic')
            dpYUp[:, x] = interp(paramRfYIndexesNew )

        #dpX
        for x in range(bl.stopX - bl.startX):
            interp = interp1d(paramRfYIndexes, dpX[:,x] , kind = 'cubic')
            dpXUp[:, x] = interp(paramRfYIndexesNew )

        #motion compensate RF data
        bl.ReadFrame()
        preRf = bl.data[bl.startY:bl.stopY - bl.stepY, bl.startX:bl.stopX].copy()
        bl.ReadFrame()
        postRf = bl.data.copy()
        
        locationY, locationX = numpy.mgrid[bl.startY:bl.stopY - bl.stepY,bl.startX:bl.stopX]  
        motionCompRf = interp2( numpy.arange(postRf.shape[1]), numpy.arange(postRf.shape[0]), postRf, locationX + dpXUp, locationY + dpYUp)

        #Now compute the cross correlation between the pre frame and the motion compensated post frame
        import cv
        template = cv.fromarray( numpy.float32(motionCompRf)  )
        image = cv.fromarray( numpy.float32( preRf)   )
        resultCv = cv.fromarray(numpy.float32( numpy.zeros( (1,1) ) ) )
        cv.MatchTemplate(template, image, resultCv, cv.CV_TM_CCORR_NORMED )
        
        resultNp = numpy.asarray(resultCv)
        rhoRfJ[skip, 0] =  float(resultNp)


        ###################
        ####Other set###### 
        ###################
        reader.SetFileName(RESULT_DIR + '/blockMatch/frame_' + str(frameNo) + '_'+ str(frameNo + skip + 1) + 'dispY.mhd')
        dpYItk = reader.Execute()
        dpY = sitk.GetArrayFromImage(dpYItk)

        reader.SetFileName(RESULT_DIR + '/blockMatch/frame_' + str(frameNo) + '_'+ str(frameNo + skip + 1) + 'dispX.mhd')
        dpXItk = reader.Execute()
        dpX = sitk.GetArrayFromImage(dpXItk)
       
        dpYUp = zeros( (bl.stopY - bl.stepY - bl.startY , bl.stopX - bl.startX) )
        dpXUp = zeros( (bl.stopY - bl.stepY - bl.startY , bl.stopX - bl.startX) )

        #Upsample displacement to same spacing as RF.
        paramRfYIndexes = numpy.arange(bl.startY, bl.stopY, bl.stepY )
        paramRfYIndexesNew = numpy.arange(bl.startY, bl.stopY - bl.stepY )

        #dpY
        for x in range(bl.stopX - bl.startX):
            interp = interp1d(paramRfYIndexes, dpY[:,x] , kind = 'cubic')
            dpYUp[:, x] = interp(paramRfYIndexesNew )

        #dpX
        for x in range(bl.stopX - bl.startX):
            interp = interp1d(paramRfYIndexes, dpX[:,x] , kind = 'cubic')
            dpXUp[:, x] = interp(paramRfYIndexesNew )

        #motion compensate RF data
        bl.ReadFrame()
        preRf = bl.data[bl.startY:bl.stopY - bl.stepY, bl.startX:bl.stopX].copy()
        bl.ReadFrame()
        postRf = bl.data.copy()
        
        locationY, locationX = numpy.mgrid[bl.startY:bl.stopY - bl.stepY,bl.startX:bl.stopX]  
        motionCompRf = interp2(numpy.arange(postRf.shape[1]), numpy.arange(postRf.shape[0]), postRf, locationX + dpXUp, locationY + dpYUp)

        #Now compute the cross correlation between the pre frame and the motion compensated post frame
        template = cv.fromarray( numpy.float32(motionCompRf)  )
        image = cv.fromarray( numpy.float32( preRf)   )
        resultCv = cv.fromarray(numpy.float32( numpy.zeros( (1,1) ) ) )
        cv.MatchTemplate(template, image, resultCv, cv.CV_TM_CCORR_NORMED )
        
        resultNp = numpy.asarray(resultCv)
        rhoRfK[0,skip] = float(resultNp) #Compute motion compensated correlation
    
    rhoRf = (rhoRfJ + rhoRfK)/2
    rhoS = zeros( (5,5) )


    for j in range(5):
        for k in range(5):
           ###Unfortunately going to have to go through several file formats to be able to compute normalized CC
           reader.SetFileName(RESULT_DIR + '/blockMatch/frame_' + str(frameNo - 1 - j) + '_' + str(frameNo) + '.mhd')
           imJ = reader.Execute()
           numpyImJ = sitk.GetArrayFromImage(imJ)

           reader.SetFileName(RESULT_DIR + '/blockMatch/frame_' + str(frameNo) + '_' + str(frameNo + k + 1) + '.mhd')
           imK = reader.Execute()
           numpyImK = sitk.GetArrayFromImage(imK)
           
           template = cv.fromarray( numpy.float32(abs(numpyImJ)  ) )
           image = cv.fromarray( numpy.float32( abs(numpyImK))  )
           resultCv = cv.fromarray(numpy.float32( numpy.zeros( (1,1) ) ) )
           cv.MatchTemplate(template, image, resultCv, cv.CV_TM_CCORR_NORMED )
        
           rhoS[j,k] = float(numpy.asarray(resultCv) )
    
    DQMarray = rhoS*rhoRf
    
    #Pick out combination with best DQM
    preSkip, postSkip = numpy.unravel_index(DQMarray.argmax(), DQMarray.shape)
    framePairs[frameNo, 0] = preSkip
    framePairs[frameNo, 1] = postSkip
    DQM[frameNo] = DQMarray.max()

    #Create composite strain image, normalizing to 1%
    reader.SetFileName(RESULT_DIR + '/blockMatch/frame_' + str(frameNo) + '_' + str(frameNo + 1 + postSkip) + '.mhd')
    strainPreItk = reader.Execute()
    strainPre = abs(sitk.GetArrayFromImage(strainPreItk))
    strainPre *= .01*strainPre.mean()

    reader.SetFileName(RESULT_DIR + '/blockMatch/frame_' + str(frameNo - 1 - preSkip) + '_' + str(frameNo) + '.mhd')
    strainPostItk = reader.Execute()
    strainPost = abs(sitk.GetArrayFromImage(strainPostItk))
    strainPost *= .01*strainPost.mean()

    strainComposite = (strainPre + strainPost)/2

    #Save composite strain as itk image
    itkIm = sitk.GetImageFromArray(strainComposite)
    itkIm.SetOrigin(strainPreItk.GetOrigin())
    itkIm.SetSpacing(strainPreItk.GetSpacing())
    writer = sitk.ImageFileWriter()

    writer.SetFileName(RESULT_DIR + '/blockMatch/frame_' + str(frameNo) + 'strainComposite.mhd' )
    writer.Execute(itkIm)

#Save winning frame triplets and DQM in arrays in case I want to pick out particular frames later
numpy.save(RESULT_DIR + '/framePairs' ,framePairs)
numpy.save(RESULT_DIR + '/DQM' ,DQM)

#Now create plots, showing B-mode image, strain with adjacent frame
#composite strain, and DQM
os.makedirs(RESULT_DIR + '/pngImages')

from scipy.signal import hilbert
reader = sitk.ImageFileReader()
from matplotlib import pyplot

for frameNo in range(5, bl.nFrames -5):

    bl.ReadFrame(frameNo)
    bMode = numpy.log10( abs( hilbert(bl.data, axis = 0)) )
    bMode -= bMode.max()
   
    reader.SetFileName(RESULT_DIR + '/blockMatch/frame_' + str(frameNo) + 'strainComposite.mhd')
    strainCompItk = reader.Execute()
    strainComp = sitk.GetArrayFromImage(strainCompItk)
    strainCompRGB = bl.CreateParametricImage(strainComp, strainCompItk.GetOrigin(), strainCompItk.GetSpacing(), inPixels = False, frameNo
    = frameNo, colormap = 'gray', vmin = 0, vmax = .01 )
    
    reader.SetFileName(RESULT_DIR + '/blockMatch/frame_' + str(frameNo) +'_'+ str(frameNo + 1 ) + '.mhd')
    strainItk = reader.Execute()
    strain = abs(sitk.GetArrayFromImage(strainItk))
    strainRGB = bl.CreateParametricImage(strain, strainItk.GetOrigin(), strainItk.GetSpacing(), inPixels = False, frameNo
    = frameNo, colormap = 'gray', vmin = 0, vmax = .01 )

    fig = pyplot.figure()
    
    plt = fig.add_subplot(2,2,1)
    plt.imshow(bMode, extent = [0, bl.fovX, bl.fovY, 0], cmap = 'gray' , vmin = -3, vmax = 0)
    
    plt = fig.add_subplot(2,2,2)
    plt.imshow(strainRGB, extent = [0, bl.fovX, bl.fovY, 0], cmap = 'gray', vmin = 0, vmax = 0.01 )
    
    plt = fig.add_subplot(2,2,3)
    plt.imshow(strainCompRGB, extent = [0, bl.fovX, bl.fovY, 0], cmap = 'gray', vmin = 0, vmax = 0.01 )
    
    plt = fig.add_subplot(2,2,4)
    plt.plot( DQM, lw = 5)
    plt.plot( frameNo, DQM[frameNo], 'ro', markersize = 13 )

    fig.savefig(RESULT_DIR + '/pngImages/frame_' + str(frameNo).zfill(3) + '.png')
