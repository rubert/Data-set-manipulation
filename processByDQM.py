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

for frameNo in range(bl.nFrames):
    for skipNo in range(4):

        bl.CreateStrainImage(preFrame = frameNo, skipFrame = skipNo, itkFileName = RESULT_DIR + '/blockMatch/frame_' + str(frameNo) +'_'
        + str(frameNo + 1 + skipNo) )


#After block-matching I compute the frame triplets with the optimum DQM
import numpy
from numpy import zeros

for frameNo in range(5, bl.nFrames -5):
    tmpDQM = zeros( (5,5) )
    
    rhoRfJ = zeros( (5,1) )
    rhoRfK = zeros( (1,5) )
    rhoRf = zeros( (5,5) )
    
    import simpleITK as sitk
    import cv
    reader = sitk.ImageFileReader()
    
    for skip in range(5):
        #load up displacement between frames
        reader.SetFileName(RESULT_DIR + 'blockMatch/frame_' + str(frameNo) + '_' str(frameNo + skip + 1) + 'dispY.mhd')
        dpYItk = reader.execute()
        dpY = sitk.GetArrayFromImage(frameJDispItk).T

        reader.SetFileName(RESULT_DIR + 'blockMatch/frame_' + str(frameNo) + '_' str(frameNo + skip + 1) + 'dispX.mhd')
        dpXItk = reader.execute()
        dpX = sitk.GetArrayFromImage(frameJDispItk).T
       
        dpYUp = zeros( (bl.stopY - bl.startY + 1, bl.stopX - bl.startX + 1) )
        dpXUp = zeros( (bl.stopY - bl.startY + 1, bl.stopX - bl.startX + 1) )

        #Upsample displacement to same spacing as RF.
        paramRfYIndexes = arange(bl.startY, bl.stopY, bl.stepY )
        paramRfYIndexesNew = arange(bl.startY, bl.stopY )

        #dpY
        for x in range(bl.stopX - bl.startX + 1):
            interp = interpolate.interp1d(paramRfYIndexes, dpY[:,x] )
            dpYUp[:, x] = interp(paramRfYIndexesNew )

        #dpX
        for x in range(bl.stopX - bl.startX + 1):
            interp = interpolate.interp1d(paramRfYIndexes, dpX[:,x] )
            dpXUp[:, x] = interp(paramRfYIndexesNew )

        #motion compensate RF data
        bl.ReadData()
        preRf = bl.data[startY:stopY, startX:stopX].copy()
        bl.ReadData()
        postRf = bl.data[startY:stopY, startX:stopX].copy()
        
        from scipy.interpolate import interp2d
        locationX, locationY = numpy.m_grid[startX:stopX, startY:stopY]  
        interp = interp2d(postRf, locationX, locationY)

        motionCompRf = interp(locationX + dpX, locationY + dpY)

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
        reader.SetFileName(RESULT_DIR + 'blockMatch/frame_' + str(frameNo) + '_' str(frameNo + skip + 1) + 'dispY.mhd')
        dpYItk = reader.execute()
        dpY = sitk.GetArrayFromImage(frameJDispItk).T

        reader.SetFileName(RESULT_DIR + 'blockMatch/frame_' + str(frameNo) + '_' str(frameNo + skip + 1) + 'dispX.mhd')
        dpXItk = reader.execute()
        dpX = sitk.GetArrayFromImage(frameJDispItk).T
       
        dpYUp = zeros( (bl.stopY - bl.startY + 1, bl.stopX - bl.startX + 1) )
        dpXUp = zeros( (bl.stopY - bl.startY + 1, bl.stopX - bl.startX + 1) )

        #Upsample displacement to same spacing as RF.
        paramRfYIndexes = arange(bl.startY, bl.stopY, bl.stepY )
        paramRfYIndexesNew = arange(bl.startY, bl.stopY )

        #dpY
        for x in range(bl.stopX - bl.startX + 1):
            interp = interpolate.interp1d(paramRfYIndexes, dpY[:,x] )
            dpYUp[:, x] = interp(paramRfYIndexesNew )

        #dpX
        for x in range(bl.stopX - bl.startX + 1):
            interp = interpolate.interp1d(paramRfYIndexes, dpX[:,x] )
            dpXUp[:, x] = interp(paramRfYIndexesNew )

        #motion compensate RF data
        bl.ReadData()
        preRf = bl.data[startY:stopY, startX:stopX].copy()
        bl.ReadData()
        postRf = bl.data[startY:stopY, startX:stopX].copy()
        
        locationX, locationY = numpy.m_grid[startX:stopX, startY:stopY]  
        interp = interp2d(postRf, locationX, locationY)

        motionCompRf = interp(locationX + dpX, locationY + dpY)

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
           imJ = simpleITK.load
           numpyImJ = simpleITK.convert

           imK = simpleITK.load
           numpyImK = simpleITK.convert
           
           rhoS[j,k] = cv.normCrossCorrelation(numpyImJ, numpyImK )
    
