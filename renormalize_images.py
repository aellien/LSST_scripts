import numpy as np
import os
from astropy.io import fits

if __name__ == '__main__':

    # Paths, lists & variables
    path_data = '/home/ellien/LSST_ICL/simulations/out1'
    path_scripts = '/home/ellien/LSST_ICL/scripts'
    path_plots = '/home/ellien/LSST_ICL/plots/out1'
    path_wavelets = '/n03data/ellien/LSST_ICL/wavelets/out1/'

    nfl = [ 'noise_00761_0000078_0.05_g_2Mpc.rebin.fits', \
            'noise_00761_0000120_0.05_g_2Mpc.rebin.fits', \
            'noise_00761_0000126_0.05_g_2Mpc.rebin.fits', \
            'noise_00761_0000174_0.05_g_2Mpc.rebin.fits', \
            'noise_00761_0000385_0.05_g_2Mpc.rebin.fits'  ]

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
        hduo.writeto( os .path.join( path_data, onf ) )
