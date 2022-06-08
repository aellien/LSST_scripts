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
import pandas as pd
from astropy.io import fits
from astropy.visualization import LinearStretch, LogStretch
from astropy.visualization import ZScaleInterval, MinMaxInterval
from astropy.visualization import ImageNormalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import SymLogNorm
from scipy.stats import sigmaclip
from skimage.measure import label, regionprops
from sklearn.utils import resample
from atom_props import *
from measure_icl_quantities import *

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def make_galaxy_catalog( oim, nf, n_levels, n_sig_gal = 50, level_gal = 3 ):

    # path, list & variables
    gal = np.zeros( (xs, ys) )

    # Read atoms
    nf = nf[:-4]
    nfl = ''.join( (nf, 'ol.*') )
    rimpath = os.path.join(path_wavelets, nfl)

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
    hduo.writeto(os.path.join( path_wavelets, ''.join( ( nf, 'results.residuals.fits') )), overwrite = True )

    hduo = fits.PrimaryHDU(gal)
    hduo.writeto(os.path.join( path_wavelets, ''.join( ( nf, 'results.gal.hardsepBCG.fits') )), overwrite = True )

    hduo = fits.PrimaryHDU(iclbcg)
    hduo.writeto(os.path.join( path_wavelets, ''.join( ( nf, 'results.iclbcg.hardsepBCG.fits') )), overwrite = True )

    hduo = fits.PrimaryHDU(rim)
    hduo.writeto(os.path.join( path_wavelets, ''.join( ( nf, 'results.rim.fits') )), overwrite = True )

    #rdc.to_fits( os.path.join( path_wavelets, ''.join( ( nf, 'results.rdc.fits') )), overwrite = True )

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

    # Read atoms
    nf = nf[:-4]
    nfl = ''.join( (nf, 'ol.*') )
    rimpath = os.path.join(path_wavelets, nfl)
    grimpath = glob.glob(rimpath)
    grimpath.sort()

    moim = np.max(oim)
    for i, it in enumerate( grimpath ):

        #print(it)
        ol = d.read_objects_from_pickle( it )
        itl = d.read_interscale_trees_from_pickle( os.path.join( path_wavelets, \
                                    ''.join( (nf, 'itl.it%03d.pkl'%(i+1)) ) ) )
        atom = np.zeros(oim.shape)
        #print(os.path.join( path_wavelets, \
        #                            ''.join( (nf, 'itl.it%03d.pkl'%(i+1)) ) ))
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

    hduo = fits.PrimaryHDU(satgal)
    hduo.writeto(os.path.join( path_wavelets, ''.join( ( nf, 'results.gal.hardsepBCG.fits') )), overwrite = True )

    hduo = fits.PrimaryHDU(iclbcg)
    hduo.writeto(os.path.join( path_wavelets, ''.join( ( nf, 'results.iclbcg.hardsepBCG.fits') )), overwrite = True )

    #rdc.to_fits( os.path.join( path_wavelets, ''.join( ( nf, 'results.rdc.fits') )), overwrite = True )

    return rdc, icl, gal, uicl, res, rim


def make_results_sizesep( oim, nfp, lvl_sep_big, n_hard_icl, rc, ricl, nf, xs, ys, n_levels ):

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

    '''
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

    # Datacube
    rdc = d.datacube( rdc_array, dctype = 'NONE', fheader = header )
    rim = np.sum( rdc_array, axis = 2 )
    res = oim - rim

    # write to fits
    hduo = fits.PrimaryHDU(res)
    hduo.writeto( nfp + 'results.residuals.fits', overwrite = True )

    hduo = fits.PrimaryHDU(icl)
    hduo.writeto( nfp + 'results.icl.sizesep.fits', overwrite = True )

    hduo = fits.PrimaryHDU(gal)
    hduo.writeto( nfp + 'results.gal.sizesep.fits', overwrite = True )

    hduo = fits.PrimaryHDU(rim)
    hduo.writeto( nfp + 'results.rim.fits', overwrite = True )

    #hduo = fits.PrimaryHDU(satgal)
    #hduo.writeto(os.path.join( path_wavelets, ''.join( ( nf, 'results.gal.hardsepBCG.fits') )), overwrite = True )

    #hduo = fits.PrimaryHDU(iclbcg)
    #hduo.writeto(os.path.join( path_wavelets, ''.join( ( nf, 'results.iclbcg.hardsepBCG.fits') )), overwrite = True )

    ##rdc.to_fits( os.path.join( path_wavelets, ''.join( ( nf, 'results.rdc.fits') )), overwrite = True )


    return rdc, icl, gal, res, rim

if __name__ == '__main__':

    # Paths, lists & variables
    path_data = '/n03data/ellien/LSST_ICL/simulations/out2/'
    path_scripts = '/home/ellien/LSST_ICL/LSST_scripts'
    path_wavelets = '/n03data/ellien/LSST_ICL/wavelets/out2/'

    dirl = ['HorizonAGN', 'Hydrangea', 'Magneticum', 'TNG-100']

    gamma = 0.5
    n_levels = 11
    lvl_sep_big = 5
    lvl_sep = 5
    n_hard_icl = 5
    pixscale = 1.6 # ''/pixel
    physscale = 1 # kpc/''
    rc = 20 # pixels, distance to center to be classified as gal
    ricl = 650 # pixels, distance to center to be classified as ICL

    flag = False
    for dir in dirl:

        image_dir = os.path.join( path_data, dir )
        image_files = glob.glob(image_dir+'/*norm.fits')

        for nf in image_files:

            print('\n%s'%nf)
            # Read files
            oimfile = os.path.join( path_data, nf )
            hdu = fits.open(oimfile)
            header = hdu[0].header
            oim = hdu[0].data

            xs, ys = oim.shape

            split = nf.split('/')
            nf = split[-1]
            nfp = os.path.join( path_wavelets, dir, 'run1', nf[:-4] )

            #rdc, icl, gal, res, rim = make_results_sizesep( oim, nfp, lvl_sep_big, n_hard_icl, rc, ricl, nf, xs, ys, n_levels )
            flux_icl_l, flux_gal_l, frac_icl_l, err_wr_icl_l, err_wr_gal_l = measure_icl_quantities_sizesep(  oim, nfp, gamma, lvl_sep_big, n_levels, n_iter = 1000, pdf = 'uniform', verbose = True )

            lowFicl, upFicl = bootstrap_error( np.array(flux_icl_l), 1000, alpha = 0.95  )
            lowFgal, upFgal = bootstrap_error( np.array(flux_gal_l), 1000, alpha = 0.95  )
            lowficl, upficl = bootstrap_error( np.array(frac_icl_l), 1000, alpha = 0.95  )

            print('SIZESEP | Flux ICL = %f +-(%f, %f), std = %f, Err_wr = %f' %(np.mean(flux_icl_l), np.mean(flux_icl_l) - lowFicl, upFicl - np.mean(flux_icl_l), np.std(flux_icl_l), np.sqrt(np.sum(np.array(err_wr_icl_l)**2))) )
            print('SIZESEP | Flux gal = %f +-(%f, %f), std = %f, Err_wr = %f' %(np.mean(flux_gal_l), np.mean(flux_gal_l) - lowFgal, upFgal - np.mean(flux_gal_l), np.std(flux_gal_l), np.sqrt(np.sum(np.array(err_wr_gal_l)**2))) )
            print('SIZESEP | Fraction ICL = %f +-(%f, %f), std = %f' %(np.mean(frac_icl_l), np.mean(frac_icl_l) - lowficl, upficl - np.mean(frac_icl_l), np.std(frac_icl_l)) )

            results_wavsep = measure_icl_quantities_wavsep( oim, nfp, gamma, lvl_sep_big, lvl_sep = lvl_sep, n_levels = n_levels, n_iter = 1000, verbose = False )

            print('WAVSEP | LVL = %d                | LVL = %d                | LVL = %d                |' %(lvl_sep - 1, lvl_sep, lvl_sep + 1))
            print('WAVSEP | Flux ICL = %1.3e | Flux ICL = %1.3e | Flux ICL = %1.3e |  ' %( results_wavsep[3], results_wavsep[0], results_wavsep[6] ) )
            print('WAVSEP | Flux gal = %1.3e | Flux gal = %1.3e | Flux gal = %1.3e |  ' %(results_wavsep[4], results_wavsep[1], results_wavsep[7] ) )
            print('WAVSEP | Fraction ICL = %1.3f | Fraction ICL = %1.3f |  Fraction ICL = %1.3f | ' %(results_wavsep[5], results_wavsep[2], results_wavsep[8] ) )

            if flag == False:
                results = pd.DataFrame( [[ dir, nf, np.mean(frac_icl_l), np.mean(frac_icl_l) - lowficl, upficl - np.mean(frac_icl_l), results_wavsep[2], results_wavsep[5], results_wavsep[8] ]], columns = [ 'dir', 'name', 'ICL fraction sizesep', 'err up', 'err low', 'ICL fraction wavsep', 'ICL fraction wavsep up', 'ICL fraction wavsep low' ])
                flag = True
            else:
                newresults = pd.DataFrame( [[ dir, nf, np.mean(frac_icl_l), np.mean(frac_icl_l) - lowficl, upficl - np.mean(frac_icl_l), results_wavsep[2], results_wavsep[5], results_wavsep[8] ]], columns = [ 'dir', 'name', 'ICL fraction sizesep', 'err up', 'err low', 'ICL fraction wavsep', 'ICL fraction wavsep up', 'ICL fraction wavsep low' ])
                results = pd.concat( [ results, newresults], ignore_index=True)

            n_coregal = 3
            #cat = make_galaxy_catalog( oim, nf, n_levels, n_sig_gal = 50, level_gal = 3 )
            cat = []
            #rdc, gal, iclbcg, res, rim = make_results_hardsepBCG( oim, nfp, lvl_sep_big, n_hard_icl, rc, ricl, nf, xs, ys, n_levels )
            #rdc, icl, gal, uicl, res, rim = make_results_hardsep( oim, nfp, lvl_sep_big, n_hard_icl, rc, ricl, nf, xs, ys, n_levels )
            #plot_dawis_results( oim, oicl, ogal, rdc, icl, gal, res, rim, path_plots )
