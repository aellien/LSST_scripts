#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Last modif: 08/2022
# Author: Amael Ellien
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Modules
import dawis as d
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import glob
import pandas as pd
from astropy.io import fits
from astropy.visualization import LinearStretch, LogStretch
from astropy.visualization import ZScaleInterval, MinMaxInterval
from astropy.visualization import ImageNormalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import sigmaclip
from skimage.measure import label, regionprops

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def convert_2D_to_1D(IMAGE, SIZE_IMAGE, L_BINS):

    CENTER     = SIZE_IMAGE / 2
    BINS       = np.linspace(0, SIZE_IMAGE / 2, L_BINS) # Bins for radial profile
    SUMBINS    = np.zeros(L_BINS - 1)

    # Create arrays of indexes
    X, Y = np.meshgrid( np.arange( - CENTER, CENTER), np.arange( - CENTER, CENTER) )

    # Create matrix of radii
    R    = np.sqrt( X**2 + Y**2 )

    # Measure radial ICL profile
    for i in range(L_BINS - 1):
        NORM       = np.size(IMAGE[ (R >= BINS[i]) & ( R < BINS[i + 1] ) ])
        SUMBINS[i] = np.sum(IMAGE[ (R >= BINS[i]) & ( R < BINS[i + 1] ) ])# / NORM

    return SUMBINS, BINS

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def measure_transition_radius(nfp, im_icl, im_bcg, n_bins, pixscale, physscale ):

    size_image = im_icl.shape[0]
    profile_icl, bins = convert_2D_to_1D(im_icl, size_image, n_bins)
    profile_bcg, bins = convert_2D_to_1D(im_bcg - im_icl, size_image, n_bins)

    plot_radial_profile(nfp, bins, pixscale, physscale, profile_icl, profile_bcg)

    try:
        bin_r_trans = np.where( profile_icl >= profile_iclbcg )[0][0]
        print(np.where(profile_icl > profile_bcg))
        r_trans_pix = bin_r_trans * size_image / 2. / n_bins
        r_trans_kpc = r_trans_pix * pixscale * physscale
    except:
        r_trans_kpc = np.inf

    return r_trans_kpc

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def plot_radial_profile(nfp, bins, pixscale, physscale, profile_icl, profile_bcg):

    fig = plt.figure()
    plt.plot(bins[:-1] * pixscale * physscale, profile_icl, label = 'icl')
    plt.plot(bins[:-1] * pixscale * physscale, profile_bcg, label = 'bcg')
    plt.yscale('log')
    plt.xscale('log')
    #plt.ylim(bottom=1E-3)
    plt.legend()
    plt.savefig( nfp + 'profile.png', format = 'png' )

    return None
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if __name__ == '__main__':

    # Paths, lists & variables
    path_data = '/home/ellien/LSST_ICL/simulations/out4/'
    path_scripts = '/home/ellien/LSST_ICL/LSST_scripts'
    path_wavelets = '/home/ellien/LSST_ICL/wavelets/out4/'

    dirl = ['TNG-100']
    schl = [ 'sizesep_200','wavsep_005']

    norm = 1E4
    n_bins = 300
    pixscale = 0.8 # ''/pixel
    physscale = 1 # kpc/''

    radii = []

    for sch in schl:

        for dir in dirl:

            image_dir = os.path.join( path_wavelets, dir )
            image_files = glob.glob(image_dir+'/run1/00099_0000004_0.05_yz_r_4Mpc_mu30.3.*.results.icl.%s.fits'%sch)

            for image in image_files:

                splt = image.split('.')
                name = splt[0]
                print('\n%s'%name)

                path_icl = splt
                path_icl[-3] = 'icl'
                path_icl = '.'.join(path_icl)
                print(path_icl)

                path_gal = splt
                path_gal[-3] = 'gal'
                path_gal = '.'.join(path_gal)
                print(path_gal)

                path_sat = splt
                path_sat[-3] = 'gal'
                path_sat[-2] = 'bcg'+sch
                path_sat = '.'.join(path_sat)
                print(path_sat)

                # Read files
                im_icl = fits.getdata(path_icl) * norm
                im_gal = fits.getdata(path_gal) * norm
                im_sat = fits.getdata(path_sat) * norm
                im_bcg = im_gal - im_sat

                size_image = im_icl.shape[0]

                #plt.figure()
                #plt.imshow(im_bcg, norm = ImageNormalize(im_bcg, interval = MinMaxInterval(), stretch = LogStretch() ), cmap = 'binary')
                #plt.show()

                r_trans_kpc = measure_transition_radius(im_icl, im_bcg, size_image, n_bins, pixscale, physscale )
                radii.append( [ name, sch, r_trans_kpc ] )

                print(r_trans_kpc)

                fig = plt.figure()
                plt.plot(bins[:-1] * pixscale * physscale, profile_icl, label = 'icl')
                plt.plot(bins[:-1] * pixscale * physscale, profile_bcg, label = 'bcg')
                plt.yscale('log')
                plt.xscale('log')
                plt.ylim(bottom=1E-3)
                plt.legend()
                plt.savefig( image[:-4] + 'profile.png', format = 'png' )
