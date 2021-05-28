import os, sys
import cv2 as cv
import numpy as np
from matplotlib import pyplot as plt
Rmin = int(sys.argv[1])
Rmax = int(sys.argv[2])

# Loads an image
dirpath = sys.argv[3]
# dirpath = '../20200806/D_8mm_N3/1/'
imgname = sys.argv[4]
# imgname = '1_026g.png'
imgfile = imgname
outdir = sys.argv[5]
files = os.listdir(dirpath)

files.sort()
for f in files:
    if imgfile in f:
        src = cv.imread(dirpath+f)
        if src is None:
            print ('Error opening image!')
        Rimg = src[:,:,2]
        Rimg = cv.medianBlur(Rimg, 3)

        rows = Rimg.shape[0]
        circles = cv.HoughCircles(Rimg, cv.HOUGH_GRADIENT, 1, 90,
                                   param1=100, param2=18,
                                   minRadius=Rmin, maxRadius=Rmax)

        if circles is not None:
            circles = np.uint16(np.around(circles))
            for i in circles[0, :]:
                center = (i[0], i[1])
                # circle center
                cv.circle(src, center, 1, (0, 100, 100), 3)
                # circle outline
                radius = i[2]
                cv.circle(src, center, radius, (255, 0, 255), 3)

        if not os.path.exists(outdir):
            os.mkdir(outdir)
        if not os.path.exists(outdir+'txt/'):
            os.mkdir(outdir+'txt/')
        if not os.path.exists(outdir+'png/'):
            os.mkdir(outdir+'png/')            
        np.savetxt(outdir+'txt/'+f[:-4]+'.txt', circles[0,:], delimiter='\t')
        cv.imwrite(outdir+'png/'+f, src)

        
    elif f.endswith(imgfile[imgfile.find('*')+1:]) and f.startswith(imgfile[:imgfile.find('*')]):
#         print("f_all:", f)             
        src = cv.imread(dirpath+f)
#         Check if image is loaded fine
        if src is None:
            print ('Error opening image!')
        Rimg = src[:,:,2]
    
#         Reduce the noise to avoid false circle detection
        Rimg = cv.medianBlur(Rimg, 3)

        ## [houghcircles]
        rows = Rimg.shape[0]
        circles = cv.HoughCircles(Rimg, cv.HOUGH_GRADIENT, 1, 90,
                                   param1=100, param2=18,
                                   minRadius=Rmin, maxRadius=Rmax)

        ## [draw]
        if circles is not None:
            circles = np.uint16(np.around(circles))
            for i in circles[0, :]:
                center = (i[0], i[1])
                # circle center
                cv.circle(src, center, 1, (0, 100, 100), 3)
                # circle outline
                radius = i[2]
                cv.circle(src, center, radius, (255, 0, 255), 3)

        if not os.path.exists(outdir):
            os.mkdir(outdir)
        if not os.path.exists(outdir+'txt/'):
            os.mkdir(outdir+'txt/')
        if not os.path.exists(outdir+'png/'):
            os.mkdir(outdir+'png/')            
        np.savetxt(outdir+'txt/'+f[:-4]+'.txt', circles[0,:], delimiter='\t')
        cv.imwrite(outdir+'png/'+f, src)
# '''