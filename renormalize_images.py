import numpy as np
import os
from astropy.io import fits
from scipy.ndimage import gaussian_filter

if __name__ == '__main__':

    # Paths, lists & variables
    path_data = '/n03data/ellien/LSST_ICL/simulations/out2/HorizonAGN'
    path_scripts = '/home/ellien/LSST_ICL/scripts'
    path_plots = '/home/ellien/LSST_ICL/plots/out2'
    path_wavelets = '/n03data/ellien/LSST_ICL/wavelets/out2/'

    nfl = [ '00761_0000001_0.05_xy_r_4Mpc_mu30.rebin.fits', \
            '00761_0000001_0.05_yz_r_4Mpc_mu30.rebin.fits' ]

    for nf in nfl:

        hdu = fits.open( os.path.join( path_data, nf ) )
        im, header = hdu[0].data, hdu[0].header
        norm = np.max(im)

        oim = np.copy(im)
        oim /= norm
        oheader = header
        oheader['NORM'] = norm

        onf = nf[:-4] + 'norm.fits'
        hduo = fits.PrimaryHDU( oim, header = oheader )
        hduo.writeto( os.path.join( path_data, onf ), overwrite = True )

        #fim = np.copy(im)
        #fim = gaussian_filter( fim, 1.0 )

        #fnf = nf[:-4] + 'gauss.fits'
        #hduo = fits.PrimaryHDU( fim, header = header )
        #hduo.writeto( os.path.join( path_data, fnf ), overwrite = True )
