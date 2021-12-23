def rebin_fits(FILENAME,XBIN=2,YBIN=2):

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
    #print(np.shape(im))

    #---------------------------------------------------------------------------
    # BINNAGE
    xedge = np.shape(im)[0]%XBIN
    yedge = np.shape(im)[1]%YBIN
    im = im[xedge:,yedge:]
    binim = np.reshape(im,(int(np.shape(im)[0]/XBIN),XBIN,int(np.shape(im)[1]/YBIN),YBIN)) 
    binim = np.mean(binim,axis=3)
    binim = np.mean(binim,axis=1)

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
