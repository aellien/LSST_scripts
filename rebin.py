import os
import sys
import glob

def rebin_fits(FILENAME,XBIN=2,YBIN=2,type='MEAN'):

    #---------------------------------------------------------------------------
    # MODULES
    from astropy.io import fits
    import pylab as plt
    import numpy as np

    #---------------------------------------------------------------------------
    # LECTURE FILE
    fitsFile = fits.open(FILENAME)
    im = np.array(fitsFile[0].data)
    #im = im[:7040,:7040]
    print('Input shape: %d x %d' %( im.shape[0], im.shape[1] ))
    print('XBIN = %d, YBIN = %d' %( XBIN, YBIN ))

    #---------------------------------------------------------------------------
    # BINNAGE
    xedge = np.shape(im)[0]%XBIN
    yedge = np.shape(im)[1]%YBIN
    im = im[xedge:,yedge:]
    binim = np.reshape(im,(int(np.shape(im)[0]/XBIN),XBIN,int(np.shape(im)[1]/YBIN),YBIN))

    if type == 'MEAN':
        binim = np.mean(binim,axis=3)
        binim = np.mean(binim,axis=1)
    elif type == 'SUM':
        binim = np.sum(binim,axis=3)
        binim = np.sum(binim,axis=1)
    print('New shape: %d x %d' %( binim.shape[0], binim.shape[1]  ))

    #---------------------------------------------------------------------------
    # NEW FILE
    oldheader = fitsFile[0].header
    newheader = oldheader
    try:
        newheader['CRPIX1'] = oldheader['CRPIX1']/float(XBIN)
        newheader['CD1_1'] = oldheader['CD1_1']*float(XBIN)
        newheader['CD1_2'] = oldheader['CD1_2']*float(XBIN)
        newheader['CRPIX2'] = oldheader['CRPIX2']/float(YBIN)
        newheader['CD2_2'] = oldheader['CD2_2']*float(YBIN)
        newheader['CD2_1'] = oldheader['CD2_1']*float(YBIN)
    except:
        pass

    hdu = fits.PrimaryHDU(binim,header=newheader)
    hdu.writeto(FILENAME[:-4]+'rebin.fits',overwrite=True)

if __name__ == '__main__':

    # Paths, lists & variables
    path_data = '/n03data/ellien/LSST_ICL/simulations/out4'
    path_scripts = '/home/ellien/LSST_ICL/scripts'

    dirl = ['HorizonAGN', 'Hydrangea', 'Magneticum', 'TNG-100']

    for dir in dirl:

        image_dir = os.path.join( path_data, dir )
        image_files = glob.glob(image_dir+'/*.fits')

        for im in image_files:

            rebin_fits( im, XBIN = 4, YBIN = 4, type = 'SUM')
