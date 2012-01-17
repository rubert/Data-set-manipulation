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
    print "Directory: " + RESULT_DIR + " Exists."
    sys.exit()

os.makedirs(RESULT_DIR)
os.makedirs(RESULT_DIR + '/blockMatch')


#Perform block-matching
from blockMatch import blockMatchClass
bl = blockMatchClass(BASE_PATH +'/'+ RFD_NAME, 'rfd', selectRoi = True)

#for frameNo in range(bl.nFrames):
for frameNo in range(6):
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

from matplotlib import pyplot
import pdb
pdb.set_trace()
DQM = zeros( bl.nFrames )

#for frameNo in range(5, bl.nFrames -5):
for frameNo in range(5, 6):

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
        motionCompRf = interp2( numpy.arange(postRf.shape[1]), numpy.arange(postRf.shape[0]), postRf, locationY + dpYUp, locationX + dpXUp)

        #Now compute the cross correlation between the pre frame and the motion compensated post frame
        import cv
        template = cv.fromarray( numpy.float32(motionCompRF)  )
        image = cv.fromarray( numpy.float32( preRf)   )
        resultCv = cv.fromarray(numpy.float32( np.zeros( (1,1) ) ) )
        cv.MatchTemplate(template, image, resultCv, cv.CV_TM_CCORR_NORMED )
        
        resultNp = np.asarray(resultCv)
        rhoRfJ[skip] =  float(resultNp)


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
        paramRfYIndexes = arange(bl.startY, bl.stopY, bl.stepY )
        paramRfYIndexesNew = arange(bl.startY, bl.stopY - bl.stepY )

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
        motionCompRf = interp2(numpy.arange(postRf.shape[1]), numpy.arange(postRf.shape[0]), postRf, locationY + dpYUp, locationX + dpXUp)

        #Now compute the cross correlation between the pre frame and the motion compensated post frame
        template = cv.fromarray( numpy.float32(motionCompRF)  )
        image = cv.fromarray( numpy.float32( preRf)   )
        resultCv = cv.fromarray(numpy.float32( np.zeros( (1,1) ) ) )
        cv.MatchTemplate(template, image, resultCv, cv.CV_TM_CCORR_NORMED )
        
        resultNp = np.asarray(resultCv)
        rhoRfK[skip] = float(resultNp) #Compute motion compensated correlation
    
    rhoRf = rhoRfJ + rhoRfK
    rhoS = zeros( (5,5) )


    for j in range(5):
        for k in range(5):
           ###Unfortunately going to have to go through several file formats to be able to compute normalized CC
           imJ = reader.SetFileName(RESULT_DIR + '/blockMatch/frame_' + str(frameNo - 1 - j) + '_' + str(frameNo) + '.mhd')
           numpyImJ = sitk.GetArrayFromImage(imJ)

           imK = reader.SetFileName(RESULT_DIR + '/blockMatch/frame_' + str(frameNo) + '_' + str(frameNo + k + 1) + '.mhd')
           numpyImK = sitk.GetArrayFromImage(imK)
           
           template = cv.fromarray( numpy.float32(abs(numpyImJ)  ) )
           image = cv.fromarray( numpy.float32( abs(numpyImK))  )
           resultCv = cv.fromarray(numpy.float32( np.zeros( (1,1) ) ) )
           cv.MatchTemplate(template, image, resultCv, cv.CV_TM_CCORR_NORMED )
        
           rhoS[j,k] = float(np.asarray(resultCv) )
    

    DQMarray = rhoS*rhoRf
    
    #Pick out combination with best DQM
    preSkip, postSkip = DQMarray.argmax()
    DQM[frameNo] = DQMarray.max()

    #Create composite strain image, normalizing to 1%
    reader.SetFileName(RESULT_DIR + '/blockMatch/frame_' + str(frameNo) + '_' + str(frameNo + 1 + postSkip) + 'dispY.mhd')
    strainPreItk = reader.Execute()
    strainPre = abs(sitk.GetArrayFromImage(strainPreItk).T)
    strainPre *= .01*strainPre.mean()

    reader.SetFileName(RESULT_DIR + '/blockMatch/frame_' + str(frameNo - 1 - preSkip) + '_' + str(frameNo) + 'dispY.mhd')
    strainPostItk = reader.Execute()
    strainPost = abs(sitk.GetArrayFromImage(strainPostItk).T)
    strainPost *= .01*strainPost.mean()

    strainComposite = (strainPre + strainPost)/2

    #Save composite strain as itk image
    itkIm = sitk.GetImageFromArray(strainComposite)
    writer = sitk.ImageFileWriter()

    writer.SetFileName(RESULT_DIR + '/blockMatch/frame_' + str(frameNo) + 'strainComposite.mhd' )
    writer.Execute(itkIm)


#Now create plots, showing B-mode image, strain with 1 skipped frame
#composite strain, and DQM
os.makedirs(RESULT_DIR + '/pngImages')

from scipy.signal import hilbert
for frameNo in range(5, bl.nFrames -5):

    bl.ReadFrame(frameNo)
    bMode = numpy.log10( abs( hilbert(bl.data, axis = 0)) )
   
    reader.SetFileName(RESULT_DIR + '/blockMatch/frame_' + str(frameNo) + 'strainComposite.mhd')
    strainCompItk = reader.Execute()
    strainComp = sitk.GetArrayFromImage(strainCompItk).T
    strainCompRGB = bl.CreateParametricImage(strainComp, strainCompItk.GetOrigin(), strainCompItk.GetSpacing(), frameNo
    = frameNo, colormap = 'gray', vmin = 0, vmax = .01 )
    

    reader.SetFileName(RESULT_DIR + '/blockMatch/frame_' + str(frameNo) +'_'+ str(frameNo + 1 ) + '.mhd')
    strainItk = reader.Execute()
    strain = abs(sitk.GetArrayFromImage(strainItk).T)
    strainRGB = bl.CreateParametricImage(strain, strainItk.GetOrigin(), strainItk.GetSpacing(), frameNo
    = frameNo, colormap = 'gray', vmin = 0, vmax = .01 )

    fig = pyplot.Figure()
    
    plt = fig.add_subplot(2,2,1)
    plt.imshow(bMode, extent = [0, bl.fovX, bl.fovY], cmap = 'gray' )
    
    plt = fig.add_subplot(2,2,2)
    plt.imshow(strainRGB, extent = [0, bl.fovX, bl.fovY], cmap = 'gray', vmin = 0, vmax = 0.01 )
    
    plt = fig.add_subplot(2,2,3)
    plt.imshow(strainCompRGB, extent = [0, bl.fovX, bl.fovY], cmap = 'gray', vmin = 0, vmax = 0.01 )
    
    plt = fig.add_subplot(2,2,4)
    plt.plot( DQM, lw = 5)
    plt.plot( frameNo, DQM[frameNo], 'r0', markersize = 13 )

    pyplot.savefig(RESULT_DIR + '/pngImages/frame_' + str(frameNo).zfill(3) + '.png')
