import dawis as d
import os
import sys

if __name__ == '__main__':

    # Paths, lists & variables
    path_data = '/home/ellien/LSST_ICL/simulations/out1'
    path_scripts = '/home/ellien/LSST_ICL/scripts'
    path_plots = '/home/ellien/LSST_ICL/plots/out2'
    path_wavelets = '/home/ellien/LSST_ICL/wavelets/out4/'

    framerate = 3

    # 'noise_00761_0000120_0.05_g_2Mpc.rebin.norm.fits', \
    nfl = [ 'noise_00761_0000078_0.05_g_2Mpc.rebin.norm.fits', \
            'noise_00761_0000126_0.05_g_2Mpc.rebin.norm.fits', \
            'noise_00761_0000174_0.05_g_2Mpc.rebin.norm.fits', \
            'noise_00761_0000385_0.05_g_2Mpc.rebin.norm.fits'  ]

    for nf in nfl:

        nfp = os.path.join( path_wavelets, nf[:-5] )
        d.make_gif( framerate, nfp )
