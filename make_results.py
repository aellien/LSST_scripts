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
import ray
import pandas as pd
from astropy.io import fits
from astropy.visualization import LinearStretch, LogStretch
from astropy.visualization import ZScaleInterval, MinMaxInterval, AsymmetricPercentileInterval
from astropy.visualization import ImageNormalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import SymLogNorm
from scipy.stats import sigmaclip
from skimage.measure import label, regionprops
from sklearn.utils import resample
from atom_props import *
from measure_icl_quantities import *
from scipy.stats import kurtosis
from measure_transition_radius import *

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def make_galaxy_catalog( oim, n_levels, n_sig_gal = 50, level_gal = 3, display = True ):

    # path, list & variables
    xs, ys = oim.shape
    gal = np.zeros( (xs, ys) )

    sigma, mean, gain = d.pg_noise_bissection( oim, max_err = 1E-3, n_sigmas = 3 )
    aim = d.anscombe_transform( oim, sigma, mean, gain )
    acdc, awdc = d.bspl_atrous( aim, n_levels )
    sdc = d.hard_threshold( awdc, n_sig_gal )
    sup = sdc.array[:,:,level_gal]
    lab = label( sup )
    reg = regionprops( lab )

    cat = []
    for r in reg:
        cat.append( [ r.centroid[1], r.centroid[0] ] )

    if display == True:

        fig, ax = plt.subplots( 1, 3 )
        ax[0].imshow( oim, norm = ImageNormalize( gal, \
                                          interval = MinMaxInterval(), \
                                          stretch = LogStretch()) )

        ax[1].imshow( gal, norm = ImageNormalize( gal, \
                                          interval = MinMaxInterval(), \
                                          stretch = LogStretch()) )
        ax[2].imshow( sup2 )

        for r in reg2:
            ax[0].plot( r.centroid[1], r.centroid[0], 'r+' )
            ax[1].plot( r.centroid[1], r.centroid[0], 'r+' )

        plt.show()

    return np.array(cat)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def make_results_hardsepBCG( oim, nfp, gamma, lvl_sep_big, n_hard_icl, rc, ricl, nf, xs, ys, n_levels ):

    # path, list & variables
    res = np.zeros( (xs, ys) )
    iclbcg = np.zeros( (xs, ys) )
    gal = np.zeros( (xs, ys) )
    rim = np.zeros( (xs, ys) )

    rdc_array = np.zeros( (xs, ys, n_levels) )

    xc = xs / 2.
    yc = ys / 2.

    # Read atoms
    ol, itl = read_image_atoms( nfp, verbose = True )
    for j, object in enumerate(ol):

        x_min, y_min, x_max, y_max = object.bbox
        lvlo = object.level
        xco = itl[j].interscale_maximum.x_max
        yco = itl[j].interscale_maximum.y_max

        flag_iclbcg = False
        flag_gal = False

        if object.level >= lvl_sep_big:

            if object.level >= n_hard_icl:

                if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) <= ricl:

                    iclbcg[ x_min : x_max, y_min : y_max ] += object.image
                    flag_iclbcg = True

            else:

                if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) <= rc:
                    iclbcg[ x_min : x_max, y_min : y_max ] += object.image
                    flag_iclbcg = True

        else:

            if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) <= rc:
                iclbcg[ x_min : x_max, y_min : y_max ] += object.image * gamma
                flag_iclbcg = True

            elif x_max - x_min >= 2**(8) or y_max - y_min >= 2**(8):

                if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) <= ricl:

                    iclbcg[ x_min : x_max, y_min : y_max ] += object.image * gamma
                    flag_iclbcg = True

            else:

                gal[ x_min : x_max, y_min : y_max ] += object.image * gamma
                flag_gal = True

        # all objects datacube
        if object.level >= lvl_sep_big:
            rdc_array[ x_min : x_max, y_min : y_max, object.level ] += object.image
        else:
            rdc_array[ x_min : x_max, y_min : y_max, object.level ] += object.image * gamma

    # Datacube
    rdc = d.datacube( rdc_array, dctype = 'NONE', fheader = header )
    rim = np.sum( rdc_array, axis = 2 )
    res = oim - rim

    # write to fits
    hduo = fits.PrimaryHDU(res)
    hduo.writeto(nfp + 'results.residuals.fits', overwrite = True )

    hduo = fits.PrimaryHDU(gal)
    hduo.writeto(nfp + 'results.gal.hardsepBCG.fits', overwrite = True )

    hduo = fits.PrimaryHDU(iclbcg)
    hduo.writeto(nfp + 'results.iclbcg.hardsepBCG.fits', overwrite = True )

    hduo = fits.PrimaryHDU(rim)
    hduo.writeto(nfp + 'results.rim.fits', overwrite = True )

    #rdc.to_fits( os.path.join( path_wavelets, ''.join( ( nf, 'results.rdc.fits') )), overwrite = True )

    return rdc, gal, iclbcg, res, rim

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def make_results_wavsep( oim, nfp, gamma, lvl_sep_big, lvl_sep, xs, ys, n_levels, plot_vignet = False ):
    '''Simple separation based on wavelet scale, given by parameter 'lvl_sep'.
    '''
    # path, list & variables
    res = np.zeros( (xs, ys) )
    icl = np.zeros( (xs, ys) )
    gal = np.zeros( (xs, ys) )
    rim = np.zeros( (xs, ys) )

    atom_gal = []
    coo_gal = []

    rdc_array = np.zeros( (xs, ys, n_levels) )

    xc = xs / 2.
    yc = ys / 2.

    # Read atoms
    ol, itl = read_image_atoms( nfp, verbose = False )
    onb = len(ol)
    filtered_onb = 0

    for j, o in enumerate(ol):

        x_min, y_min, x_max, y_max = o.bbox
        k = kurtosis(o.image.flatten(), fisher=True)
        if k < 0:
            filtered_onb += 1
            continue

        lvlo = o.level
        xco = itl[j].interscale_maximum.x_max
        yco = itl[j].interscale_maximum.y_max

        if o.level >= lvl_sep_big:

            if o.level >= lvl_sep:
                icl[ x_min : x_max, y_min : y_max ] += o.image
            else:
                gal[ x_min : x_max, y_min : y_max ] += o.image

        else:

            if o.level >= lvl_sep:
                icl[ x_min : x_max, y_min : y_max ] += o.image * gamma
            else:
                gal[ x_min : x_max, y_min : y_max ] += o.image * gamma

        if o.level >= lvl_sep_big:
            rdc_array[ x_min : x_max, y_min : y_max, o.level ] += o.image
        else:
            rdc_array[ x_min : x_max, y_min : y_max, o.level ] += o.image * gamma

    print('Kurtosis filtered: %d/%d'%(filtered_onb,onb))

    # Datacube
    rdc = d.datacube( rdc_array, dctype = 'NONE' )
    rim = np.sum( rdc_array, axis = 2 )
    res = oim - rim

    # write to fits
    hduo = fits.PrimaryHDU(res)
    hduo.writeto( nfp + 'results.residuals.fits', overwrite = True )

    hduo = fits.PrimaryHDU(icl)
    hduo.writeto( nfp + 'results.icl.wavsep_%03d.fits'%lvl_sep, overwrite = True )

    hduo = fits.PrimaryHDU(gal)
    hduo.writeto( nfp + 'results.gal.wavsep_%03d.fits'%lvl_sep, overwrite = True )

    hduo = fits.PrimaryHDU(rim)
    hduo.writeto( nfp + 'results.rim.fits', overwrite = True )

    if plot_vignet == True:
        interval = AsymmetricPercentileInterval(5, 99.5) # meilleur rendu que MinMax or ZScale pour images reconstruites
        fig, ax = plt.subplots(1, 2)
        poim = ax[0].imshow(gal, norm = ImageNormalize( gal, interval = interval, stretch = LogStretch()), cmap = 'binary')
        divider = make_axes_locatable(ax[0])
        cax = divider.append_axes("top", size="5%", pad=0.05)
        caxre = fig.colorbar( poim, cax = cax, \
                                    orientation = 'horizontal', \
                                    format = '%2.1f',\
                                    pad = 0,\
                                    shrink = 1.0,\
                                    ticklocation = 'top' )
        poim = ax[1].imshow(icl, norm = ImageNormalize( icl, interval = interval, stretch = LogStretch()), cmap = 'binary')
        divider = make_axes_locatable(ax[1])
        cax = divider.append_axes("top", size="5%", pad=0.05)
        caxre = fig.colorbar( poim, cax = cax, \
                                    orientation = 'horizontal', \
                                    format = '%2.1f',\
                                    pad = 0,\
                                    shrink = 1.0,\
                                    ticklocation = 'top' )
        plt.tight_layout()
        plt.savefig( nfp + 'results.wavsep_%03d.png'%lvl_sep, format = 'png' )
        print('Write vignet to' + nfp + 'results.wavsep_%03d.png'%(lvl_sep))
        plt.close('all')

    return icl, gal

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def make_results_sizesep( oim, nfp, gamma, lvl_sep_big, size_sep, size_sep_pix, xs, ys, n_levels, plot_vignet = False ):

    # Paths, list & variables
    atom_icl = []
    atom_gal = []
    xs, ys = oim.shape
    im_icl = np.zeros((xs, ys))
    im_gal = np.zeros((xs, ys))

    # Read atoms and interscale trees
    ol, itl = read_image_atoms( nfp, verbose = False )
    for j, o in enumerate(ol):

        x_min, y_min, x_max, y_max = o.bbox
        k = kurtosis(o.image.flatten(), fisher=True)
        if k < 0:
            continue

        sx = x_max - x_min
        sy = y_max - y_min

        if o.level >= lvl_sep_big:
            if (sx >= size_sep_pix) or (sy >= size_sep_pix):
                atom_icl.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
                im_icl[ x_min : x_max, y_min : y_max ] += o.image
            else:
                atom_gal.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
                im_gal[ x_min : x_max, y_min : y_max ] += o.image
        else:
            if (sx >= size_sep_pix) or (sy >= size_sep_pix):
                atom_icl.append( [ x_max - x_min, y_max - y_min, np.sum(o.image) * gamma, o.level ] )
                im_icl[ x_min : x_max, y_min : y_max ] += o.image * gamma
            else:
                atom_gal.append( [ x_max - x_min, y_max - y_min, np.sum(o.image) * gamma, o.level ] )
                im_gal[ x_min : x_max, y_min : y_max ] += o.image * gamma

    atom_icl = np.array(atom_icl)
    atom_gal = np.array(atom_gal)

    hduo = fits.PrimaryHDU(im_icl)
    hduo.writeto( nfp + 'results.icl.sizesep_%03d.fits'%size_sep, overwrite = True )

    hduo = fits.PrimaryHDU(im_gal)
    hduo.writeto( nfp + 'results.gal.sizesep_%03d.fits'%size_sep, overwrite = True )

    if plot_vignet == True:

        fig, ax = plt.subplots(1, 2)
        ax[0].imshow(im_gal, norm = ImageNormalize(im_gal, interval = MinMaxInterval(), stretch = LogStretch() ), cmap = 'binary')
        ax[1].imshow(im_icl, norm = ImageNormalize(im_icl, interval = MinMaxInterval(), stretch = LogStretch() ), cmap = 'binary')
        plt.tight_layout()
        plt.savefig( nfp + 'results.sizesep_%03d.png'%size_sep, format = 'png' )
        plt.close('all')

    return im_icl, im_gal

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def make_results_sbt( oim, nfp, gamma, lvl_sep_big, sbt, norm, pixscale, xs, ys, n_levels, plot_vignet = False ):

    # Paths, list & variables
    atom_icl = []
    atom_gal = []
    xs, ys = oim.shape
    im_tot = np.zeros((xs, ys))

    # Read atoms and interscale trees
    ol, itl = read_image_atoms( nfp, verbose = False )

    # Atom wavelet scale separation
    for j, o in enumerate(ol):
        x_min, y_min, x_max, y_max = o.bbox
        k = kurtosis(o.image.flatten(), fisher=True)
        if k < 0:
            continue

        if o.level >= lvl_sep_big:
            im_tot[ x_min : x_max, y_min : y_max ] += o.image
        else:
            im_tot[ x_min : x_max, y_min : y_max ] += o.image * gamma

    #im_tot *= norm # renormalize for SB
    sbt_flux = 10**(-0.4 * sbt) * pixscale**2 / norm # renormalize for SB
    sb_lim = 10**(-0.4 * 30.3) * pixscale**2 / norm # renormalize for SB

    print(sbt_flux * norm)

    im_icl = np.copy(im_tot)
    im_icl[im_icl >= sbt_flux] = 0.
    im_icl[im_icl < sb_lim] = 0.

    im_gal = np.copy(im_tot)
    im_gal[im_gal <= sbt_flux] = 0.


    '''
    im_tot[im_tot < 10E-30] = 10E-30 # get rid of nul pixels
    im_tot_sb = - 2.5 * np.log10(im_tot / pixscale**2 )

    im_icl = np.copy(im_tot_sb)
    im_icl[im_icl <= sbt] = 0.
    im_icl[im_icl > 30.3] = 0.
    im_gal = np.copy(im_tot_sb)
    im_gal[im_gal > sbt] = 0.
    '''

    hduo = fits.PrimaryHDU(im_icl)
    hduo.writeto( nfp + 'results.icl.sbt_%2.1f.fits'%sbt, overwrite = True )

    hduo = fits.PrimaryHDU(im_gal)
    hduo.writeto( nfp + 'results.gal.sbt_%2.1f.fits'%sbt, overwrite = True )

    if plot_vignet == True:
        fig, ax = plt.subplots(1, 2)
        im_gal[im_gal == 0. ] = sbt
        im_icl[im_icl == 0. ] = sbt
        ax[0].imshow(im_gal, norm = ImageNormalize(im_gal, interval = MinMaxInterval(), vmin = im_gal.min(), vmax = sbt, stretch = LinearStretch() ), cmap = 'binary_r')
        ax[1].imshow(im_icl, norm = ImageNormalize(im_icl, interval = MinMaxInterval(), vmin = sbt, vmax = 35, stretch = LinearStretch() ), cmap ='binary_r')
        plt.tight_layout()
        plt.savefig( nfp + 'results.sbt_%2.1f.png'%sbt, format = 'png' )
        plt.close('all')

    return im_icl, im_gal

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def make_results_spawavsep( oim, nfp, gamma, lvl_sep_big, cat_gal, rc, n_sig_gal, lvl_sep, xs, ys, n_levels, plot_vignet = False  ):
    '''
    Wavelet + spatial separation. ICL and ICL+BCG.
    '''
    # path, list & variables
    icl = np.zeros( (xs, ys) )
    gal = np.zeros( (xs, ys) )

    atom_gal = []
    coo_gal = []

    rdc_array = np.zeros( (xs, ys, n_levels) )

    xc = xs / 2.
    yc = ys / 2.

    # Read atoms
    ol, itl = read_image_atoms( nfp, verbose = True )

    for j, o in enumerate(ol):

        x_min, y_min, x_max, y_max = o.bbox
        lvlo = o.level
        xco = itl[j].interscale_maximum.x_max
        yco = itl[j].interscale_maximum.y_max

        if o.level >= lvl_sep_big:

            if o.level < lvl_sep:

                flag_gal = False
                for pos in cat_gal:

                    if np.sqrt( ( xco - pos[1] )**2 + ( yco - pos[0] )**2 ) <= rc:
                        flag_gal = True
                        break

                if flag_gal == True:
                    if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) <= rc:
                        icl[ x_min : x_max, y_min : y_max ] += o.image
                    else:
                        gal[ x_min : x_max, y_min : y_max ] += o.image

                if flag_gal == False:
                    icl[ x_min : x_max, y_min : y_max ] += o.image

            else:
                icl[ x_min : x_max, y_min : y_max ] += o.image

        else:

            if o.level < lvl_sep:

                flag_gal = False
                for pos in cat_gal:

                    if np.sqrt( ( xco - pos[1] )**2 + ( yco - pos[0] )**2 ) <= rc:
                        flag_gal = True
                        break

                if flag_gal == True:
                    if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) <= rc:
                        icl[ x_min : x_max, y_min : y_max ] += o.image * gamma
                    else:
                        gal[ x_min : x_max, y_min : y_max ] += o.image * gamma

                if flag_gal == False:
                    icl[ x_min : x_max, y_min : y_max ] += o.image * gamma

            else:
                icl[ x_min : x_max, y_min : y_max ] += o.image * gamma

    hduo = fits.PrimaryHDU(icl)
    hduo.writeto( nfp + 'results.iclbcg.spawavsep_%03d.fits'%lvl_sep, overwrite = True )

    hduo = fits.PrimaryHDU(gal)
    hduo.writeto( nfp + 'results.gal.spawavsep_%03d.fits'%lvl_sep, overwrite = True )

    if plot_vignet == True:

        fig, ax = plt.subplots(1, 2)
        ax[0].imshow(gal, norm = ImageNormalize(gal, vmax = np.max(gal) / 10., interval = MinMaxInterval(), stretch = LogStretch() ), cmap = 'binary')
        ax[1].imshow(icl, norm = ImageNormalize(icl, vmax = np.max(icl) / 10., interval = MinMaxInterval(), stretch = LogStretch() ), cmap = 'binary')
        plt.tight_layout()
        plt.savefig( nfp + 'results.spawavsep_%03d.png'%lvl_sep, format = 'png' )
        plt.close('all')

    return icl, gal

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def make_results_bcgwavsep( oim, nfp, gamma, lvl_sep_big, rc, lvl_sep, xs, ys, n_levels, plot_vignet = False  ):
    '''
    Wavelet + spatial separation. ICL and ICL+BCG.
    '''
    # path, list & variables
    icl = np.zeros( (xs, ys) )
    gal = np.zeros( (xs, ys) )

    atom_gal = []
    coo_gal = []

    rdc_array = np.zeros( (xs, ys, n_levels) )

    xc = xs / 2.
    yc = ys / 2.

    # Read atoms
    ol, itl = read_image_atoms( nfp, verbose = False )

    for j, o in enumerate(ol):

        x_min, y_min, x_max, y_max = o.bbox
        k = kurtosis(o.image.flatten(), fisher=True)
        if k < 0:
            continue

        lvlo = o.level
        xco = itl[j].interscale_maximum.x_max
        yco = itl[j].interscale_maximum.y_max

        if o.level >= lvl_sep_big:

            if o.level >= lvl_sep:
                icl[ x_min : x_max, y_min : y_max ] += o.image
            else:
                if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) <= rc:
                    icl[ x_min : x_max, y_min : y_max ] += o.image
                else:
                    gal[ x_min : x_max, y_min : y_max ] += o.image

        else:

            if o.level >= lvl_sep:
                icl[ x_min : x_max, y_min : y_max ] += o.image * gamma
            else:
                if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) <= rc:
                    icl[ x_min : x_max, y_min : y_max ] += o.image * gamma
                else:
                    gal[ x_min : x_max, y_min : y_max ] += o.image * gamma

    hduo = fits.PrimaryHDU(icl)
    hduo.writeto( nfp + 'results.iclbcg.bcgwavsep_%03d.fits'%lvl_sep, overwrite = True )

    hduo = fits.PrimaryHDU(gal)
    hduo.writeto( nfp + 'results.gal.bcgwavsep_%03d.fits'%lvl_sep, overwrite = True )

    if plot_vignet == True:

        fig, ax = plt.subplots(1, 2)
        ax[0].imshow(gal, norm = ImageNormalize(gal, vmax = np.max(gal) / 10., interval = MinMaxInterval(), stretch = LogStretch() ), cmap = 'binary')
        ax[1].imshow(icl, norm = ImageNormalize(icl, vmax = np.max(icl) / 10., interval = MinMaxInterval(), stretch = LogStretch() ), cmap = 'binary')
        plt.tight_layout()
        plt.savefig( nfp + 'results.bcgwavsep_%03d.png'%lvl_sep, format = 'png' )
        plt.close('all')

    return icl, gal

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def make_results_bcgsizesep( oim, nfp, gamma, lvl_sep_big, rc, size_sep, size_sep_pix, xs, ys, n_levels, plot_vignet = False ):

    # Paths, list & variables
    atom_icl = []
    atom_gal = []
    xs, ys = oim.shape
    im_icl = np.zeros((xs, ys))
    im_gal = np.zeros((xs, ys))
    xc = xs / 2.
    yc = ys / 2.

    # Read atoms and interscale trees
    ol, itl = read_image_atoms( nfp, verbose = False )
    for j, o in enumerate(ol):

        x_min, y_min, x_max, y_max = o.bbox
        k = kurtosis(o.image.flatten(), fisher=True)
        if k < 0:
            continue

        sx = x_max - x_min
        sy = y_max - y_min
        xco = itl[j].interscale_maximum.x_max
        yco = itl[j].interscale_maximum.y_max

        if o.level >= lvl_sep_big:
            if (sx >= size_sep_pix) or (sy >= size_sep_pix):
                atom_icl.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
                im_icl[ x_min : x_max, y_min : y_max ] += o.image
            else:
                if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) <= rc: #BCG
                    im_icl[ x_min : x_max, y_min : y_max ] += o.image
                else:
                    atom_gal.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
                    im_gal[ x_min : x_max, y_min : y_max ] += o.image
        else:
            if (sx >= size_sep_pix) or (sy >= size_sep_pix):
                atom_icl.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
                im_icl[ x_min : x_max, y_min : y_max ] += o.image * gamma
            else:
                if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) <= rc: #BCG
                    im_icl[ x_min : x_max, y_min : y_max ] += o.image * gamma
                else:
                    atom_gal.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
                    im_gal[ x_min : x_max, y_min : y_max ] += o.image * gamma

    atom_icl = np.array(atom_icl)
    atom_gal = np.array(atom_gal)

    hduo = fits.PrimaryHDU(im_icl)
    hduo.writeto( nfp + 'results.iclbcg.bcgsizesep_%03d.fits'%size_sep, overwrite = True )

    hduo = fits.PrimaryHDU(im_gal)
    hduo.writeto( nfp + 'results.gal.bcgsizesep_%03d.fits'%size_sep, overwrite = True )

    if plot_vignet == True:

        fig, ax = plt.subplots(1, 2)
        ax[0].imshow(im_gal, norm = ImageNormalize(im_gal, vmax = np.max(im_gal) / 10., interval = MinMaxInterval(), stretch = LogStretch() ), cmap = 'binary')
        ax[1].imshow(im_icl, norm = ImageNormalize(im_icl, vmax = np.max(im_icl) / 10., interval = MinMaxInterval(), stretch = LogStretch() ), cmap = 'binary')
        plt.tight_layout()
        plt.savefig( nfp + 'results.bcgsizesep_%03d.png'%size_sep, format = 'png' )
        plt.close('all')

    return im_icl, im_gal

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def make_results_sizesep2( oim, nfp, gamma, lvl_sep_big, rc, size_sep_bcg, size_sep_icl, size_sep_bcg_pix, size_sep_icl_pix, xs, ys, n_levels, plot_vignet = False ):
    '''
    02/2023: results for ICL only from sizesep (make_results_sizesep) are not good when compared to simulation values. Too much BCG flux in ICL.
    --> new methodology: create BCG+ICL map, and then remove atoms smaller than size_sep_bcg_pix.
    '''
    # Paths, list & variables
    atom_icl = []
    atom_gal = []
    xs, ys = oim.shape
    im_icl = np.zeros((xs, ys))
    im_gal = np.zeros((xs, ys))
    xc = xs / 2.
    yc = ys / 2.

    # Read atoms and interscale trees
    ol, itl = read_image_atoms( nfp, verbose = False )
    for j, o in enumerate(ol):

        x_min, y_min, x_max, y_max = o.bbox
        sx = x_max - x_min
        sy = y_max - y_min
        xco = itl[j].interscale_maximum.x_max
        yco = itl[j].interscale_maximum.y_max

        if o.level >= lvl_sep_big:
            if (sx >= size_sep_icl_pix) or (sy >= size_sep_icl_pix): # isolate ICL+BCG from satellite galaxies
                if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) <= rc: #if BCG
                    if (sx >= size_sep_bcg_pix) or (sy >= size_sep_bcg_pix): #rm BCG from ICL+BCG
                        im_icl[ x_min : x_max, y_min : y_max ] += o.image
                        atom_icl.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
                    else:
                        atom_gal.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
                        im_gal[ x_min : x_max, y_min : y_max ] += o.image
                else:
                    im_icl[ x_min : x_max, y_min : y_max ] += o.image
                    atom_icl.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
            else:
                atom_gal.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
                im_gal[ x_min : x_max, y_min : y_max ] += o.image
        else:
            if (sx >= size_sep_icl_pix) or (sy >= size_sep_icl_pix): # isolate ICL+BCG from satellite galaxies
                if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) <= rc: #if BCG
                    if (sx >= 2*size_sep_bcg_pix) or (sy >= 2*size_sep_bcg_pix): #rm BCG from ICL+BCG
                        im_icl[ x_min : x_max, y_min : y_max ] += o.image * gamma
                        atom_icl.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
                    else:
                        atom_gal.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
                        im_gal[ x_min : x_max, y_min : y_max ] += o.image * gamma
                else:
                    im_icl[ x_min : x_max, y_min : y_max ] += o.image * gamma
                    atom_icl.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
            else:
                atom_gal.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
                im_gal[ x_min : x_max, y_min : y_max ] += o.image * gamma

    atom_icl = np.array(atom_icl)
    atom_gal = np.array(atom_gal)

    hduo = fits.PrimaryHDU(im_icl)
    hduo.writeto( nfp + 'results.icl.sizesep2_%03d.fits'%size_sep_icl, overwrite = True )

    hduo = fits.PrimaryHDU(im_gal)
    hduo.writeto( nfp + 'results.gal.sizesep2_%03d.fits'%size_sep_icl, overwrite = True )

    if plot_vignet == True:

        fig, ax = plt.subplots(1, 2)
        ax[0].imshow(im_gal, norm = ImageNormalize(im_gal, vmax = np.max(im_gal) / 10., interval = MinMaxInterval(), stretch = LogStretch() ), cmap = 'binary')
        ax[1].imshow(im_icl, norm = ImageNormalize(im_icl, vmax = np.max(im_icl) / 10., interval = MinMaxInterval(), stretch = LogStretch() ), cmap = 'binary')
        plt.tight_layout()
        plt.savefig( nfp + 'results.bcgsizesep_%03d.png'%size_sep, format = 'png' )
        plt.close('all')

    return im_icl, im_gal

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def trash():

    '''
        # path, list & variables
        res = np.zeros( (xs, ys) )
        icl = np.zeros( (xs, ys) )
        uicl = np.zeros( (xs, ys) )
        gal = np.zeros( (xs, ys) )
        rim = np.zeros( (xs, ys) )
        unclass = np.zeros( (xs, ys) )

        lres = []
        licl = []
        luicl = []
        lgal = []
        lrim = []
        lunclass = []

        atom_gal = []
        coo_gal = []

        rdc_array = np.zeros( (xs, ys, n_levels) )

        xc = xs / 2.
        yc = ys / 2.

        ol, itl = read_image_atoms( nfp, verbose = True )
        df_sizes = average_size_atom( ol, n_levels )
        print(df_sizes)

        lvl_sx = np.argmax( df_sizes['dsx_n'][2:] ) + 2
        lvl_sy = np.argmax( df_sizes['dsy_n'][2:] ) + 2

        sx = df_sizes['<sx>'][lvl_sx]
        if sx == 0.:
            sx = df_sizes['<sx>'][lvl_sx +1]
        sy = df_sizes['<sy>'][lvl_sy]
        if sy == 0.:
            sy = df_sizes['<sx>'][lvl_sy +1]
        print(lvl_sx, sx)
        print(lvl_sy, sy)

        for j, object in enumerate(ol):

            x_min, y_min, x_max, y_max = object.bbox
            lvlo = object.level
            xco = itl[j].interscale_maximum.x_max
            yco = itl[j].interscale_maximum.y_max

            flag_icl = False
            flag_gal = False
            flag_uicl = False

            # ICL or gal based on size
            if x_max - x_min >= sx or y_max - y_min >= sy:

                if object.level >= lvl_sep_big:
                    icl[ x_min : x_max, y_min : y_max ] += object.image
                else:
                    icl[ x_min : x_max, y_min : y_max ] += object.image * gamma

            else:

                if object.level >= lvl_sep_big:
                    gal[ x_min : x_max, y_min : y_max ] += object.image
                else:
                    gal[ x_min : x_max, y_min : y_max ] += object.image * gamma
                atom_gal.append(object)
                coo_gal.append([xco, yco])

            # All object dc
            if object.level >= lvl_sep_big:
                rdc_array[ x_min : x_max, y_min : y_max, object.level ] += object.image
            else:
                rdc_array[ x_min : x_max, y_min : y_max, object.level ] += object.image * gamma


        bcg = np.copy(gal)

        for i, ogal in enumerate(atom_gal):

            x_min, y_min, x_max, y_max = ogal.bbox
            xco, yco = coo_gal[i]
            lvlo = ogal.level

            if lvlo <= 4:

                if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) > rc * 2**(lvlo-1):

                    if x_max - x_min <= 300 or y_max - y_min <= 300:

                        bcg[ x_min : x_max, y_min : y_max ] -= ogal.image * gamma

                    elif np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) > ricl:

                        bcg[ x_min : x_max, y_min : y_max ] -= ogal.image * gamma

        satgal = np.copy(gal)
        satgal -= bcg
        iclbcg = np.copy(icl)
        iclbcg += bcg

    '''
    '''
                if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) <= ricl:

                    icl[ x_min : x_max, y_min : y_max ] += object.image
                    flag_icl = True

                    if object.level >= n_hard_icl + 1:
                        uicl[ x_min : x_max, y_min : y_max ] += object.image
            else:
                gal[ x_min : x_max, y_min : y_max ] += object.image
                atom_gal.append(object)
                coo_gal.append([xco, yco])
                ## galaxies
                #for pos in cat:
                #    if np.sqrt( ( xco - pos[1] )**2 + ( yco - pos[0] )**2 ) <= rc:
                #        gal1[ x_min : x_max, y_min : y_max ] += object.image * gamma
                #        flag_gal = True
                #        break
                #if flag_gal == False:
                #    icl[ x_min : x_max, y_min : y_max ] += object.image * gamma
        else:
            #if x_max - x_min >= 2**(n_hard_icl) or y_max - y_min >= 2**(n_hard_icl):
            #    if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) <= ricl:
            #        icl[ x_min : x_max, y_min : y_max ] += object.image * gamma
            #        flag_icl = True
            #else:
            gal[ x_min : x_max, y_min : y_max ] += object.image * gamma
            atom_gal.append(object)
            coo_gal.append([xco, yco])
            flag_gal = True
        # all objects datacube
        if object.level >= lvl_sep_big:
            rdc_array[ x_min : x_max, y_min : y_max, object.level ] += object.image
        else:
            rdc_array[ x_min : x_max, y_min : y_max, object.level ] += object.image * gamma


    bcg = np.copy(gal)

    for i, ogal in enumerate(atom_gal):

        x_min, y_min, x_max, y_max = ogal.bbox
        xco, yco = coo_gal[i]
        lvlo = ogal.level

        if lvlo <= 4:

            if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) > rc * 2**(lvlo-1):

                if x_max - x_min <= 300 or y_max - y_min <= 300:

                    bcg[ x_min : x_max, y_min : y_max ] -= ogal.image * gamma

                elif np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) > ricl:

                    bcg[ x_min : x_max, y_min : y_max ] -= ogal.image * gamma

    satgal = np.copy(gal)
    satgal -= bcg
    iclbcg = np.copy(icl)
    iclbcg += bcg
    '''
    #iclbcg = np.copy(icl)
    #satgal = np.copy(gal)

    #sizes = average_size_atom(atom_gal, n_levels)

    #for i, ogal in enumerate(atom_gal):

    #    x_min, y_min, x_max, y_max = ogal.bbox
    #    xco, yco = coo_gal[i]
    #    lvlo = ogal.level

    #    if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) <= rc * 2**(lvlo-1):

    #        if ogal.level >= lvl_sep_big:
    #            iclbcg[ x_min : x_max, y_min : y_max ] += ogal.image
    #            satgal[ x_min : x_max, y_min : y_max ] -= ogal.image
    #        else:
    #            iclbcg[ x_min : x_max, y_min : y_max ] += ogal.image * gamma
    #            satgal[ x_min : x_max, y_min : y_max ] -= ogal.image * gamma

    #    elif x_max - x_min >= 2**(n_hard_icl) or y_max - y_min >= 2**(n_hard_icl):

    #        if ogal.level >= lvl_sep_big:
    #            iclbcg[ x_min : x_max, y_min : y_max ] += ogal.image
    #            satgal[ x_min : x_max, y_min : y_max ] -= ogal.image
    #        else:
    #            iclbcg[ x_min : x_max, y_min : y_max ] += ogal.image * gamma
    #            satgal[ x_min : x_max, y_min : y_max ] -= ogal.image * gamma
    return None

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
@ray.remote
def make_results_cluster( oim, nfp, dir, nf, xs, ys, gamma, n_levels, lvl_sep_big, lvl_sep_l, size_sep_icl_l, size_sep_bcg_l, sbt_l, err_size, pixscale, physscale, rc, rc_pix, n_sig_gal, ricl, r_lsst):
    '''
    Runs all classification schemes for a single cluster. Performed by a single ray worker.
    '''

    # WAVSEP
    for lvl_sep in lvl_sep_l:

        icl, gal = make_results_wavsep( oim, nfp, gamma, lvl_sep_big, lvl_sep, xs, ys, n_levels, plot_vignet = True )
        results_wavsep = measure_icl_quantities_wavsep( oim, nfp, gamma = gamma, lvl_sep_big = lvl_sep_big, lvl_sep = lvl_sep, n_levels = n_levels, r_lsst = r_lsst, verbose = False )

    # SIZESEP
    for size_sep in size_sep_icl_l:

        size_sep_pix = size_sep * 2. / pixscale * physscale
        icl, gal = make_results_sizesep( oim, nfp, gamma, lvl_sep_big, size_sep, size_sep_pix, xs, ys, n_levels, plot_vignet = True )
        results_sizesep = measure_icl_quantities_sizesep( oim, nfp, gamma = gamma, size_sep = size_sep_pix, err_size = err_size, lvl_sep_big = lvl_sep_big, n_levels = n_levels, r_lsst = r_lsst, verbose = False )

    # SIZESEP 2 TEST
    #for size_sep in size_sep_l:

    #    size_sep_icl = size_sep
    #    size_sep_bcg = size_sep + 100. # premier test radius kpc

    #    size_sep_icl_pix = size_sep_icl * 2. / pixscale * physscale # separation diameter in pixel
    #    size_sep_bcg_pix = size_sep_bcg * 2. / pixscale * physscale # separation diameter in pixel

    #    icl, gal = make_results_sizesep2( oim, nfp, gamma, lvl_sep_big, ricl, size_sep_bcg, size_sep_icl, size_sep_bcg_pix, size_sep_icl_pix, xs, ys, n_levels, plot_vignet = False )
    #    print('sizesep2', np.sum(icl)/(np.sum(icl)+np.sum(gal) ))
    #    #results_sizesep = measure_icl_quantities_sizesep( oim, nfp, gamma = gamma, size_sep = size_sep_pix, err_size = err_size, lvl_sep_big = lvl_sep_big, n_levels = n_levels, r_lsst = r_lsst, verbose = False )
    #    #print('SIZESEP | %12s%9d | %12s%9d | %12s%9d |' %('SIZE_LOW = ', size_sep * ( 1 - err_size ) , 'SIZE = ', size_sep, 'SIZE_UP = ', size_sep * ( 1 + err_size ) ))
    #    #print('SIZESEP | %12s%1.3e | %12s%1.3e | %12s%1.3e |' %( 'Flux ICL = ', results_sizesep[6], 'Flux ICL = ', results_sizesep[0], 'Flux ICL = ', results_sizesep[3] ) )
    #    #print('SIZESEP | %12s%1.3e | %12s%1.3e | %12s%1.3e |  ' %('Flux gal = ', results_sizesep[7], 'Flux gal = ', results_sizesep[1], 'Flux gal = ', results_sizesep[4] ) )
    #    #print('SIZESEP | %12s%1.3e | %12s%1.3e | %12s%1.3e | ' %('fICL = ', results_sizesep[8], 'fICL = ', results_sizesep[2], 'fICL = ', results_sizesep[5] ) )

    # SBT
    for sbt in sbt_l:

        norm = header['NORM']
        icl, gal = make_results_sbt(oim, nfp, gamma, lvl_sep_big, sbt, norm, pixscale, xs, ys, n_levels, plot_vignet = True)
        results_sbt = measure_icl_quantities_sbt( oim, nfp, gamma = gamma, pixscale = pixscale, lvl_sep_big = lvl_sep_big, sbt = sbt, norm = norm, n_levels = n_levels, r_lsst = r_lsst, verbose = False  )

    # BCGWAVSEP
    for lvl_sep in lvl_sep_l:

        icl, gal = make_results_bcgwavsep( oim, nfp, gamma, lvl_sep_big, rc_pix, lvl_sep, xs, ys, n_levels, plot_vignet = True )
        results_bcgwavsep = measure_icl_quantities_bcgwavsep( oim, nfp, gamma, lvl_sep_big, rc_pix, lvl_sep, xs, ys, n_levels, r_lsst = r_lsst, verbose = False )

    # BCGSIZESEP
    for size_sep in size_sep_bcg_l:
        size_sep_pix = size_sep * 2. / pixscale * physscale
        icl, gal = make_results_bcgsizesep( oim, nfp, gamma, lvl_sep_big, rc_pix, size_sep, size_sep_pix, xs, ys, n_levels, plot_vignet = True )
        results_bcgsizesep = measure_icl_quantities_bcgsizesep( oim, nfp, gamma, lvl_sep_big, rc_pix, size_sep_pix, xs, ys, n_levels, r_lsst, verbose = False )

    results = pd.DataFrame( [[ dir, nf, results_wavsep[0] * norm, results_wavsep[1] * norm, results_wavsep[2], results_sizesep[0] * norm, results_sizesep[1] * norm, results_sizesep[2], \
                                results_sbt[0], results_sbt[1], results_sbt[2], results_bcgwavsep[0] * norm, results_bcgwavsep[1] * norm, results_bcgwavsep[2], \
                                results_bcgsizesep[0] * norm, results_bcgsizesep[1] * norm, results_bcgsizesep[2] ]], \
                columns = [ 'dir', 'name', 'FICL_wavsep', 'Fgal_wavsep', 'fICL_wavsep', 'FICL_sizesep', 'Fgal_sizesep', 'fICL_sizesep', \
                                'FICL_sbt', 'Fgal_sbt', 'fICL_sbt', 'FICL_bcgwavsep', 'Fgal_bcgwavsep', 'fICL_bcgwavsep', 'FICL_bcgsizesep', 'Fgal_bcgsizesep', 'fICL_bcgsizesep' ])

    return results


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if __name__ == '__main__':

    # Paths, lists & variables
    path_data = '/n03data/ellien/LSST_ICL/simulations/out4/'
    path_scripts = '/home/ellien/LSST_ICL/LSST_scripts'
    path_wavelets = '/n03data/ellien/LSST_ICL/wavelets/out4/'

    dirl = ['HorizonAGN', 'Hydrangea', 'Magneticum', 'TNG-100']
    #dirl = ['TNG-100']

    # Classification parameters
    gamma = 0.8
    n_levels = 11
    lvl_sep_big = 6
    lvl_sep_l = [ 5 ]
    size_sep_icl_l = [ 200 ] # separation radius gal/icl kpc
    size_sep_bcg_l = [ 60 ] # separation radius bcg/icl kpc
    sbt_l = [ 26. ]# [  26, 26.5, 27, 27.5, 28. ]
    err_size = 0.2
    pixscale = 0.8 # ''/pixel
    physscale = 1 # kpc/''
    rc = 30 # kpc, distance to center to be classified as gal
    rc_pix = rc  / physscale / pixscale # pixels
    n_sig_gal = 50
    ricl = 1000 # pixels, distance to center to be classified as ICL
    r_lsst = 1000 / physscale / pixscale # pixels, ICL measured within this radius (LSST choice)

    # Ray parameters
    ray_refs = []
    ray_outputs = []
    n_cpus = 48
    ray.init(num_cpus = n_cpus)
    id_gamma = ray.put(gamma)
    id_n_levels = ray.put(n_levels)
    id_lvl_sep_big = ray.put(lvl_sep_big)
    id_lvl_sep_l = ray.put(lvl_sep_l)
    id_size_sep_icl_l = ray.put(size_sep_icl_l)
    id_size_sep_bcg_l = ray.put(size_sep_bcg_l)
    id_sbt_l = ray.put(sbt_l)
    id_err_size = ray.put(err_size)
    id_pixscale = ray.put(pixscale)
    id_physscale = ray.put(physscale)
    id_rc = ray.put(rc)
    id_rc_pix = ray.put(rc_pix)
    id_n_sig_gal =  ray.put(n_sig_gal)
    id_ricl = ray.put(ricl)
    id_r_lsst = ray.put(r_lsst)

    # iterate over cluster images
    flag_data = False
    for dir in dirl:

        image_dir = os.path.join( path_data, dir )
        image_files = glob.glob(image_dir+'/*norm.fits')

        names = []
        for image in image_files:
            names.append(image.split('.')[0])
        names = np.unique(names)

        #image_files = [ '/n03data/ellien/LSST_ICL/simulations/out4/TNG-100/00099_0000004_0.05_xz_r_4Mpc_mu30.3.rebin.norm.fits', \
        #                '/n03data/ellien/LSST_ICL/simulations/out4/TNG-100/00099_0000004_0.05_xy_r_4Mpc_mu30.3.rebin.norm.fits', \
        #                '/n03data/ellien/LSST_ICL/simulations/out4/TNG-100/00099_0000004_0.05_yz_r_4Mpc_mu30.3.rebin.norm.fits' ]

        for nf in image_files:

            print('\n%s'%nf)
            # Read files
            oimfile = os.path.join( path_data, dir, nf )
            hdu = fits.open(oimfile)
            header = hdu[0].header
            oim = hdu[0].data

            xs, ys = oim.shape

            split = nf.split('/')
            nf = split[-1]
            nfp = os.path.join( path_wavelets, dir, 'run1', nf[:-4] )

            id_nf = ray.put(nf)
            id_dir = ray.put(dir)

            ray_refs.append( make_results_cluster.remote( oim,\
                                                          nfp, \
                                                          id_dir,\
                                                          id_nf,\
                                                          xs,\
                                                          ys,\
                                                          id_gamma,\
                                                          id_n_levels,\
                                                          id_lvl_sep_big,\
                                                          id_lvl_sep_l,\
                                                          id_size_sep_icl_l,\
                                                          id_size_sep_bcg_l,\
                                                          id_sbt_l,\
                                                          id_err_size,\
                                                          id_pixscale,\
                                                          id_physscale,\
                                                          id_rc,\
                                                          id_rc_pix,\
                                                          id_n_sig_gal,\
                                                          id_ricl,\
                                                          id_r_lsst ) )
    for ref in ray_refs:
        ray_outputs.append(ray.get(ref))

    ray.shutdown()

    results = ray_outputs[0]
    for output in ray_outputs[1:]:
        results = pd.concat( [ results, output], ignore_index = True )

    results.to_excel('/home/ellien/LSST_ICL/catalogs/dawis_lsst_results_%03d_%03d_%03d.xlsx'%(lvl_sep_l[0], size_sep_l[0], sbt_l[0]), sheet_name = 'amael')
