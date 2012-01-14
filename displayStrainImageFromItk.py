#!/usr/bin/python

from sys import argv,exit
import os
if len(argv) < 2:
    print 'Error, usage is:'
    print os.path.basename(argv[0]) + ': strainFile + [B mode file] + [maxStrain]'
    exit()

strainFile = argv[1]
bModeFile = None
maxStrain = None

if len(argv) == 3:
    #'figure out whether I have a B-mode file or a max strain'

    try:
        maxStrain =float(argv[2])
    except:
        bModeFile = argv[2]
        maxStrain = None

if len(argv) == 4:
    #'figure out which argument is B-mode file and which is max strain'

    try:
        maxStrain =float(argv[2])
        bModeFile = argv[3]
    except:
        maxStrain = float(argv[3])
        bModeFile = argv[2]


import SimpleITK as sitk
from matplotlib import pyplot
reader = sitk.ImageFileReader()

if bModeFile:
        
        if bModeFile == '.mhd':
            #read in B-mode from ITK file
            reader.SetFileName(bModeFile)
            bModeItk = reader.Execute()
            bMode = sitk.GetArrayFromImage(bModeItk).T
            bModeSpacing = bModeItk.GetSpacing()
            bModeOrigin = bModeItk.GetOrigin()

        else:
            #Read in RF data, let the origin be 0,0
            from rfData import rfClass
            rf = rfClass(bModeFile, 'rfd')
            rf.ReadFrame(0)
            temp = rf.data

            from scipy.signal import hilbert
            from numpy import log10
            bMode = log10(abs(hilbert(temp, axis = 0)))
            bMode = bMode - bMode.max()
            bModeOrigin = [0., 0.]
            bModeSpacing = [rf.deltaY, rf.deltaX]

        bMode[bMode < -3] = -3

        reader.SetFileName(strainFile)
        paramImItk = reader.Execute()
        paramImage = sitk.GetArrayFromImage(paramImItk).T
        paramSpacing = paramImItk.GetSpacing()
        paramOrigin = paramImItk.GetOrigin()
       
        if maxStrain:
            paramImage[paramImage > maxStrain] = maxStrain

        print "Mean value in parametric image is: " + str(paramImage.mean() )

        spacing = [ int( round(paramSpacing[0]/bModeSpacing[0]) ), int( round(paramSpacing[1]/bModeSpacing[1]) ) ]
        origin = [int( round(paramOrigin[0]/bModeSpacing[0]) ),int( round(paramOrigin[1]/bModeSpacing[1]) ) ]

        points,lines = bMode.shape

        #work out size of region in B-mode image
        bModeSizeY = spacing[0]*(paramImage.shape[0] - 1) + 1
        bModeSizeX = spacing[1]*(paramImage.shape[1] - 1) + 1
        
        from numpy import arange,zeros
        from scipy import interpolate

        paramImageUpInY = zeros( (bModeSizeY, paramImage.shape[1] ) )
        paramImageUp = zeros( (bModeSizeY, bModeSizeX) )

        #Upsample strain image to same spacing as B-mode image.
        paramRfYIndexes = arange(origin[0], origin[0] + spacing[0]*paramImage.shape[0], spacing[0] )
        paramRfYIndexesNew = arange(origin[0], origin[0] + spacing[0]*(paramImage.shape[0]-1) + 1 )
        paramRfXIndexes = arange(origin[1], origin[1] + spacing[1]*paramImage.shape[1], spacing[1] )
        paramRfXIndexesNew = arange(origin[1], origin[1] + spacing[1]*(paramImage.shape[1]-1) + 1 )

        #use old image to interpolate
        for x in range(paramImage.shape[1]):
            interp = interpolate.interp1d(paramRfYIndexes, paramImage[:,x], kind = 'nearest' )
            paramImageUpInY[:, x] = interp(paramRfYIndexesNew )


        for y in range(paramImageUp.shape[0]):
            interp = interpolate.interp1d(paramRfXIndexes, paramImageUpInY[y,:] )
            paramImageUp[y, :] = interp( paramRfXIndexesNew )


        '''Convert array containing param values to RGBALpha array'''
        from matplotlib import cm

        palette = cm.ScalarMappable()
        palette.set_cmap('gray')
        tempIm = palette.to_rgba(paramImageUp)

        palette = cm.ScalarMappable()
        palette.set_cmap('gray')
        paramBmode = palette.to_rgba(bMode)
        paramBmode[origin[0]:origin[0] + tempIm.shape[0],origin[1]:origin[1] + tempIm.shape[1], :] = tempIm

        fig =pyplot.figure()
        ax = fig.add_subplot(121)
        pyplot.imshow( paramBmode, extent = [0, bModeSpacing[1]*lines, bModeSpacing[0]*points, 0] , vmin
        =paramImage.min(), vmax = paramImage.max() )

        ax = fig.add_subplot(122)
        pyplot.imshow( bMode, extent = [0, bModeSpacing[1]*lines, bModeSpacing[0]*points, 0] , cmap = 'gray')
        pyplot.show()

else:
    #No B-mode image for reference, just display the strain
        reader.SetFileName(strainFile)
        paramImItk = reader.Execute()
        paramImage = sitk.GetArrayFromImage(paramImItk).T
        paramSpacing = paramImItk.GetSpacing()
        paramOrigin = paramImItk.GetOrigin()
        '''Convert array containing param values to RGBALpha array'''

        fig =pyplot.figure()
        if maxStrain: 
            pyplot.imshow( paramImage, extent = [0, paramSpacing[1]*paramImage.shape[1],
            paramSpacing[0]*paramImage.shape[0], 0] , cmap = 'gray', vmax = maxStrain)
        else:
            pyplot.imshow( paramImage, extent = [0, paramSpacing[1]*paramImage.shape[1],
            paramSpacing[0]*paramImage.shape[0], 0] , cmap = 'gray')
        pyplot.colorbar()
            
        pyplot.show()
