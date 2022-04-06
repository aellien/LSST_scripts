#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Last modif: 08/2021
# Author: Amael Ellien
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Modules
import dawis as d
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import glob
from astropy.io import fits
from astropy.visualization import LinearStretch, LogStretch
from astropy.visualization import ZScaleInterval, MinMaxInterval
from astropy.visualization import ImageNormalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import SymLogNorm
from scipy.stats import sigmaclip
from skimage.measure import label, regionprops

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def make_galaxy_catalog( oim, nf, n_levels, n_sig_gal = 50, level_gal = 3 ):

    # path, list & variables
    gal = np.zeros( (xs, ys) )

    # Read atoms
    nf = nf[:-4]
    nfl = ''.join( (nf, 'ol.*') )
    rimpath = os.path.join(path_wavelets, nfl)

    '''
    moim = np.max(oim)

    # make image
    for it in glob(rimpath):

        ol = d.read_objects_from_pickle( it )
        atom = np.zeros(oim.shape)

        for object in ol:

            x_min, y_min, x_max, y_max = object.bbox
            xco = x_min + ( x_max - x_min ) / 2
            yco = y_min + ( y_max - y_min ) / 2

            if np.max(object.image) > moim:
                continue

            if object.level <= n_coregal:
                gal[ x_min : x_max, y_min : y_max ] += object.image * gamma

    # get catalog

    sup = np.zeros( gal.shape )
    sup[ np.where( gal > 1E-11 ) ] = 1.

    lab = label( sup )
    reg = regionprops( lab )

    plt.ion()
    fig, ax = plt.subplots( 1, 3 )
    ax[0].imshow( oim, norm = ImageNormalize( gal, \
                                      interval = MinMaxInterval(), \
                                      stretch = LogStretch()) )

    ax[1].imshow( gal, norm = ImageNormalize( gal, \
                                      interval = MinMaxInterval(), \
                                      stretch = LogStretch()) )
    ax[2].imshow( sup )

    for r in reg:
        ax[0].plot( r.centroid[1], r.centroid[0], 'r+' )
        ax[1].plot( r.centroid[1], r.centroid[0], 'r+' )
    '''

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

    #fig, ax = plt.subplots( 1, 3 )
    #ax[0].imshow( oim, norm = ImageNormalize( gal, \
    #                                  interval = MinMaxInterval(), \
    #                                  stretch = LogStretch()) )

    #ax[1].imshow( gal, norm = ImageNormalize( gal, \
    #                                  interval = MinMaxInterval(), \
    #                                  stretch = LogStretch()) )
    #ax[2].imshow( sup2 )

    #for r in reg2:
    #    ax[0].plot( r.centroid[1], r.centroid[0], 'r+' )
    #    ax[1].plot( r.centroid[1], r.centroid[0], 'r+' )

    return np.array(cat)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def make_results_hardsepBCG( oim, path_wavelets, lvl_sep_big, n_hard_icl, rc, ricl, nf, xs, ys, n_levels ):

    # path, list & variables
    res = np.zeros( (xs, ys) )
    iclbcg = np.zeros( (xs, ys) )
    gal = np.zeros( (xs, ys) )
    rim = np.zeros( (xs, ys) )

    rdc_array = np.zeros( (xs, ys, n_levels) )

    xc = xs / 2.
    yc = ys / 2.

    # Read atoms
    nf = nf[:-4]
    nfl = ''.join( (nf, 'ol.*') )
    rimpath = os.path.join(path_wavelets, nfl)
    grimpath = glob.glob(rimpath)
    grimpath.sort()

    moim = np.max(oim)
    for i, it in enumerate( grimpath ):

        print(it)
        ol = d.read_objects_from_pickle( it )
        itl = d.read_interscale_trees_from_pickle( os.path.join( path_wavelets, \
                                    ''.join( (nf, 'itl.it%03d.pkl'%(i+1)) ) ) )
        atom = np.zeros(oim.shape)
        print(os.path.join( path_wavelets, \
                                    ''.join( (nf, 'itl.it%03d.pkl'%(i+1)) ) ))
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

                elif x_max - x_min >= 2**(n_hard_icl) or y_max - y_min >= 2**(n_hard_icl):

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
    hduo.writeto(os.path.join( path_wavelets, ''.join( ( nf, 'results.residuals.fits') )), overwrite = True )

    hduo = fits.PrimaryHDU(gal)
    hduo.writeto(os.path.join( path_wavelets, ''.join( ( nf, 'results.gal.hardsepBCG.fits') )), overwrite = True )

    hduo = fits.PrimaryHDU(iclbcg)
    hduo.writeto(os.path.join( path_wavelets, ''.join( ( nf, 'results.iclbcg.hardsepBCG.fits') )), overwrite = True )

    hduo = fits.PrimaryHDU(rim)
    hduo.writeto(os.path.join( path_wavelets, ''.join( ( nf, 'results.rim.fits') )), overwrite = True )

    rdc.to_fits( os.path.join( path_wavelets, ''.join( ( nf, 'results.rdc.fits') )), overwrite = True )

    return rdc, gal, iclbcg, res, rim


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def make_results_hardsep( oim, path_wavelets, lvl_sep_big, n_hard_icl, rc, ricl, nf, xs, ys, n_levels ):

    # path, list & variables
    res = np.zeros( (xs, ys) )
    icl = np.zeros( (xs, ys) )
    uicl = np.zeros( (xs, ys) )
    gal = np.zeros( (xs, ys) )
    rim = np.zeros( (xs, ys) )
    unclass = np.zeros( (xs, ys) )

    rdc_array = np.zeros( (xs, ys, n_levels) )

    xc = xs / 2.
    yc = ys / 2.

    # Read atoms
    nf = nf[:-4]
    nfl = ''.join( (nf, 'ol.*') )
    rimpath = os.path.join(path_wavelets, nfl)
    grimpath = glob.glob(rimpath)
    grimpath.sort()

    moim = np.max(oim)
    for i, it in enumerate( grimpath ):

        print(it)
        ol = d.read_objects_from_pickle( it )
        itl = d.read_interscale_trees_from_pickle( os.path.join( path_wavelets, \
                                    ''.join( (nf, 'itl.it%03d.pkl'%(i+1)) ) ) )
        atom = np.zeros(oim.shape)
        print(os.path.join( path_wavelets, \
                                    ''.join( (nf, 'itl.it%03d.pkl'%(i+1)) ) ))
        for j, object in enumerate(ol):

            x_min, y_min, x_max, y_max = object.bbox
            lvlo = object.level
            xco = itl[j].interscale_maximum.x_max
            yco = itl[j].interscale_maximum.y_max

            flag_icl = False
            flag_gal = False
            flag_uicl = False

            if object.level >= lvl_sep_big:

                if object.level >= n_hard_icl:

                    if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) <= ricl:

                        icl[ x_min : x_max, y_min : y_max ] += object.image
                        flag_icl = True

                        if object.level >= n_hard_icl + 1:
                            uicl[ x_min : x_max, y_min : y_max ] += object.image

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
    hduo.writeto(os.path.join( path_wavelets, ''.join( ( nf, 'results.residuals.fits') )), overwrite = True )

    hduo = fits.PrimaryHDU(icl)
    hduo.writeto(os.path.join( path_wavelets, ''.join( ( nf, 'results.icl.hardsep.fits') )), overwrite = True )

    hduo = fits.PrimaryHDU(gal)
    hduo.writeto(os.path.join( path_wavelets, ''.join( ( nf, 'results.gal.hardsep.fits') )), overwrite = True )

    hduo = fits.PrimaryHDU(uicl)
    hduo.writeto(os.path.join( path_wavelets, ''.join( ( nf, 'results.uicl.hardsep.fits') )), overwrite = True )

    hduo = fits.PrimaryHDU(rim)
    hduo.writeto(os.path.join( path_wavelets, ''.join( ( nf, 'results.rim.fits') )), overwrite = True )

    rdc.to_fits( os.path.join( path_wavelets, ''.join( ( nf, 'results.rdc.fits') )), overwrite = True )

    return rdc, icl, gal, uicl, res, rim

if __name__ == '__main__':

    # Paths, lists & variables
    path_data = '/home/ellien/LSST_ICL/simulations/out2/'
    path_scripts = '/home/ellien/LSST_ICL/scripts'
    path_wavelets = '/home/ellien/LSST_ICL/wavelets/out2/'

    dirl = ['Hydrangea', 'Magneticum', 'TNG-100'] #'HorizonAGN'

    gamma = 0.5
    n_levels = 11
    lvl_sep_big = 6
    n_hard_icl = 6
    pixscale = 1.6 # ''/pixel
    physscale = 1 # kpc/''
    rc = 20 # pixels, distance to center to be classified as gal
    ricl = 5000 # pixels, distance to center to be classified as ICL

    for dir in dirl:

        image_dir = os.path.join( path_data, dir )
        image_files = glob.glob(image_dir+'/*.fits')

        for nf in image_files:

            # Read files
            oimfile = os.path.join( path_data, nf )
            hdu = fits.open(oimfile)
            header = hdu[0].header
            oim = hdu[0].data

            xs, ys = oim.shape

            split = nf.split('/')
            nf = split[-1]
            nfp = os.path.join( path_wavelets, dir, 'run1' )

            n_coregal = 3
            #cat = make_galaxy_catalog( oim, nf, n_levels, n_sig_gal = 50, level_gal = 3 )
            cat = []
            rdc, gal, iclbcg, res, rim = make_results_hardsepBCG( oim, nfp, lvl_sep_big, n_hard_icl, rc, ricl, nf, xs, ys, n_levels )
            rdc, icl, gal, uicl, res, rim = make_results_hardsep( oim, nfp, lvl_sep_big, n_hard_icl, rc, ricl, nf, xs, ys, n_levels )
            #plot_dawis_results( oim, oicl, ogal, rdc, icl, gal, res, rim, path_plots )
