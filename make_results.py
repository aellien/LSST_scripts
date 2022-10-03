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
from matplotlib.colors import SymLogNorm
from scipy.stats import sigmaclip
from skimage.measure import label, regionprops
from sklearn.utils import resample
from atom_props import *
from measure_icl_quantities import *

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def make_galaxy_catalog( oim, nf, n_levels, n_sig_gal = 50, level_gal = 3, dislay = True ):

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
def make_results_hardsepBCG( oim, nfp, lvl_sep_big, n_hard_icl, rc, ricl, nf, xs, ys, n_levels ):

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
def make_results_wavsep( oim, nfp, lvl_sep_big, lvl_sep, xs, ys, n_levels, plot_vignet = False ):
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
    ol, itl = read_image_atoms( nfp, verbose = True )

    for j, object in enumerate(ol):

        x_min, y_min, x_max, y_max = object.bbox
        lvlo = object.level
        xco = itl[j].interscale_maximum.x_max
        yco = itl[j].interscale_maximum.y_max

        if object.level >= lvl_sep_big:

            if object.level >= lvl_sep:
                icl[ x_min : x_max, y_min : y_max ] += object.image
            else:
                gal[ x_min : x_max, y_min : y_max ] += object.image

        else:

            if object.level >= lvl_sep:
                icl[ x_min : x_max, y_min : y_max ] += object.image * gamma
            else:
                gal[ x_min : x_max, y_min : y_max ] += object.image * gamma

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
    hduo.writeto( nfp + 'results.residuals.fits', overwrite = True )

    hduo = fits.PrimaryHDU(icl)
    hduo.writeto( nfp + 'results.icl.wavsep_%d.fits'%lvl_sep, overwrite = True )

    hduo = fits.PrimaryHDU(gal)
    hduo.writeto( nfp + 'results.gal.wavsep_%d.fits'%lvl_sep, overwrite = True )

    hduo = fits.PrimaryHDU(rim)
    hduo.writeto( nfp + 'results.rim.fits', overwrite = True )

    if plot_vignet == True:

       fig, ax = plt.subplots(1, 2)
       ax[0].imshow(gal, norm = ImageNormalize(gal, interval = MinMaxInterval(), stretch = LogStretch() ), cmap = 'binary')
       ax[1].imshow(icl, norm = ImageNormalize(icl, interval = MinMaxInterval(), stretch = LogStretch() ), cmap = 'binary')
       plt.tight_layout()
       plt.savefig( nfp + 'results.wavsep_%d.fits'%lvl_sep, format = 'pdf' )
        plt.close('all')
    return None

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def make_results_sizesep( oim, nfp, lvl_sep_big, size_sep, xs, ys, n_levels, plot_vignet = False ):

    # Paths, list & variables
    atom_icl = []
    atom_gal = []
    xs, ys = oim.shape
    im_icl = np.zeros((xs, ys))
    im_gal = np.zeros((xs, ys))

    # Read atoms and interscale trees
    ol, itl = read_image_atoms( nfp, verbose = True )
    for j, o in enumerate(ol):

        x_min, y_min, x_max, y_max = o.bbox
        sx = x_max - x_min
        sy = y_max - y_min

        if o.level >= lvl_sep_big:
            if (sx >= size_sep) or (sy >= size_sep):
                atom_icl.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
                im_icl[ x_min : x_max, y_min : y_max ] += o.image
            else:
                atom_gal.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
                im_gal[ x_min : x_max, y_min : y_max ] += o.image
        else:
            if (sx >= size_sep) or (sy >= size_sep):
                atom_icl.append( [ x_max - x_min, y_max - y_min, np.sum(o.image) * gamma, o.level ] )
                im_icl[ x_min : x_max, y_min : y_max ] += o.image * gamma
            else:
                atom_gal.append( [ x_max - x_min, y_max - y_min, np.sum(o.image) * gamma, o.level ] )
                im_gal[ x_min : x_max, y_min : y_max ] += o.image * gamma

    atom_icl = np.array(atom_icl)
    atom_gal = np.array(atom_gal)

    hduo = fits.PrimaryHDU(im_icl)
    hduo.writeto( nfp + 'results.icl.sizesep_%d.fits'%size_sep, overwrite = True )

    hduo = fits.PrimaryHDU(im_gal)
    hduo.writeto( nfp + 'results.gal.sizesep_%d.fits'%size_sep, overwrite = True )

    if plot_vignet == True:

       fig, ax = plt.subplots(1, 2)
       ax[0].imshow(gal, norm = ImageNormalize(gal, interval = MinMaxInterval(), stretch = LogStretch() ), cmap = 'binary')
       ax[1].imshow(icl, norm = ImageNormalize(icl, interval = MinMaxInterval(), stretch = LogStretch() ), cmap = 'binary')
       plt.tight_layout()
       plt.savefig( nfp + 'results.sizesep_%d.fits'%size_sep, format = 'pdf' )
        plt.close('all')
    return

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def make_results_sbt( oim, nfp, lvl_sep_big, sbt, norm, xs, ys, n_levels, plot_vignet = False ):

    # Paths, list & variables
    atom_icl = []
    atom_gal = []
    xs, ys = oim.shape
    im_tot = np.zeros((xs, ys))

    # Read atoms and interscale trees
    ol, itl = read_image_atoms( nfp, verbose = True )

    # Atom wavelet scale separation
    for j, o in enumerate(ol):
        x_min, y_min, x_max, y_max = o.bbox
        if o.level >= lvl_sep_big:
            im_tot[ x_min : x_max, y_min : y_max ] += o.image
        else:
            im_tot[ x_min : x_max, y_min : y_max ] += o.image * gamma

    im_tot *= norm # renormalize for SB
    im_tot[im_tot < 10E-30] = 10E-30 # get rid of nul pixels
    im_tot_sb = - 2.5 * np.log10(im_tot / pixscale**2 )

    im_icl = np.copy(im_tot_sb)
    im_icl[im_icl <= sbt] = 0.
    im_gal = np.copy(im_tot_sb)
    im_gal[im_gal > sbt] = 0.

    hduo = fits.PrimaryHDU(im_icl)
    hduo.writeto( nfp + 'results.icl.sbt_%2.1f.fits'%sbt, overwrite = True )

    hduo = fits.PrimaryHDU(im_gal)
    hduo.writeto( nfp + 'results.gal.sbt_%2.1f.fits'%sbt, overwrite = True )

    if plot_vignet == True:
        fig, ax = plt.subplots(1, 2)
        gal[gal == 0. ] = sbt
        icl[icl == 0. ] = sbt
        ax[0].imshow(gal, norm = ImageNormalize(gal, interval = MinMaxInterval(), vmin = gal.min(), vmax = sbt, stretch = LinearStretch() ), cmap = 'binary_r')
        ax[1].imshow(icl, norm = ImageNormalize(icl, interval = MinMaxInterval(), vmin = sbt, vmax = 35, stretch = LinearStretch() ), cmap ='binary_r')
        plt.tight_layout()
        plt.savefig( nfp + 'results.sbt_%2.1f.fits'%sbt, format = 'pdf' )
        plt.close('all')
    return None

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
if __name__ == '__main__':

    # Paths, lists & variables
    path_data = '/n03data/ellien/LSST_ICL/simulations/out4/'
    path_scripts = '/home/ellien/LSST_ICL/LSST_scripts'
    path_wavelets = '/n03data/ellien/LSST_ICL/wavelets/out4/'

    dirl = ['HorizonAGN', 'Hydrangea', 'Magneticum', 'TNG-100']

    gamma = 0.8
    n_levels = 11
    lvl_sep_big = 6
    lvl_sep_l = [ 3, 4, 5, 6, 7 ]
    size_sep_l = [ 20, 40, 60, 80, 100 ] # kpc
    sbt_l = [ 26, 26.5, 27, 27.5, 28. ]
    err_size = 0.2
    pixscale = 0.8 # ''/pixel
    physscale = 1 # kpc/''
    rc = 40 # pixels, distance to center to be classified as gal
    ricl = 1000 # pixels, distance to center to be classified as ICL
    r_lsst = 1000 * physscale / pixscale # pixels, ICL measured within this radius (LSST choice)

    flag = False
    for dir in dirl:

        image_dir = os.path.join( path_data, dir )
        image_files = glob.glob(image_dir+'/*norm.fits')

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

            #-------------------------------------------------------------------
            # WAVSEP
            for lvl_sep in lvl_sep_l:

                #rdc, icl, gal, res, rim = make_results_sizesep( oim, nfp, lvl_sep_big, lvl_sep, rc, ricl, nf, xs, ys, n_levels )
                #flux_icl_l, flux_gal_l, frac_icl_l, err_wr_icl_l, err_wr_gal_l = measure_icl_quantities_sizesep(  oim, nfp, gamma, lvl_sep_big, n_levels, n_iter = 1000, pdf = 'uniform', verbose = True )

                #lowFicl, upFicl = bootstrap_error( np.array(flux_icl_l), 1000, alpha = 0.95  )
                #lowFgal, upFgal = bootstrap_error( np.array(flux_gal_l), 1000, alpha = 0.95  )
                #lowficl, upficl = bootstrap_error( np.array(frac_icl_l), 1000, alpha = 0.95  )

                #print('SIZESEP | Flux ICL = %f +-(%f, %f), std = %f, Err_wr = %f' %(np.mean(flux_icl_l), np.mean(flux_icl_l) - lowFicl, upFicl - np.mean(flux_icl_l), np.std(flux_icl_l), np.sqrt(np.sum(np.array(err_wr_icl_l)**2))) )
                #print('SIZESEP | Flux gal = %f +-(%f, %f), std = %f, Err_wr = %f' %(np.mean(flux_gal_l), np.mean(flux_gal_l) - lowFgal, upFgal - np.mean(flux_gal_l), np.std(flux_gal_l), np.sqrt(np.sum(np.array(err_wr_gal_l)**2))) )
                #print('SIZESEP | Fraction ICL = %f +-(%f, %f), std = %f' %(np.mean(frac_icl_l), np.mean(frac_icl_l) - lowficl, upficl - np.mean(frac_icl_l), np.std(frac_icl_l)) )

                make_results_wavsep( oim, nfp, lvl_sep_big, lvl_sep, xs, ys, n_levels, plot_vignet = True )
                results_wavsep = measure_icl_quantities_wavsep( oim, nfp, gamma = gamma, lvl_sep_big = lvl_sep_big, lvl_sep = lvl_sep, n_levels = n_levels, r_lsst = r_lsst, verbose = False )

                print('WAVSEP | %12s%9d | %12s%9d | %12s%9d |' %('LVL = ', lvl_sep - 1, 'LVL = ',lvl_sep, 'LVL = ',lvl_sep + 1))
                print('WAVSEP | %12s%1.3e | %12s%1.3e | %12s%1.3e |' %( 'Flux ICL = ', results_wavsep[3], 'Flux ICL = ', results_wavsep[0], 'Flux ICL = ', results_wavsep[6] ) )
                print('WAVSEP | %12s%1.3e | %12s%1.3e | %12s%1.3e |  ' %('Flux gal = ', results_wavsep[4], 'Flux gal = ', results_wavsep[1], 'Flux gal = ', results_wavsep[7] ) )
                print('WAVSEP | %12s%1.3e | %12s%1.3e | %12s%1.3e | ' %('fICL = ', results_wavsep[5], 'fICL = ', results_wavsep[2], 'fICL = ', results_wavsep[8] ) )

                #if flag == False:
                #    results = pd.DataFrame( [[ dir, nf, np.mean(frac_icl_l), np.mean(frac_icl_l) - lowficl, upficl - np.mean(frac_icl_l), \
                #results_wavsep[2], results_wavsep[5], results_wavsep[8] ]], columns = [ 'dir', 'name', 'ICL fraction sizesep', 'err up', 'err low', 'ICL fraction wavsep', 'ICL fraction wavsep up', 'ICL fraction wavsep low' ])
                #    flag = True
                #else:
                #    newresults = pd.DataFrame( [[ dir, nf, np.mean(frac_icl_l), np.mean(frac_icl_l) - lowficl, upficl - np.mean(frac_icl_l), results_wavsep[2], results_wavsep[5], results_wavsep[8] ]], columns = [ 'dir', 'name',\
                # 'ICL fraction sizesep', 'err up', 'err low', 'ICL fraction wavsep', 'ICL fraction wavsep up', 'ICL fraction wavsep low' ])
                #    results = pd.concat( [ results, newresults], ignore_index=True)

            #-------------------------------------------------------------------
            # SIZESEP
            for size_sep in size_sep_l:
                size_sep_pix = size_sep * 2. / pixscale * physscale
                make_results_sizesep( oim, nfp, lvl_sep_big, size_sep_pix, xs, ys, n_levels, plot_vignet = True )
                results_wavsep = measure_icl_quantities_sizesep( oim, nfp, gamma = gamma, size_sep = size_sep_pix, err_size = err_size, lvl_sep_big = lvl_sep_big, n_levels = n_levels, r_lsst = r_lsst, verbose = False )
                print('SIZESEP | %12s%9d | %12s%9d | %12s%9d |' %('SIZE_LOW = ', size_sep * ( 1 - err_size ) , 'SIZE = ', size_sep, 'SIZE_UP = ', size_sep * ( 1 + err_size ) ))
                print('SIZESEP | %12s%1.3e | %12s%1.3e | %12s%1.3e |' %( 'Flux ICL = ', results_wavsep[6], 'Flux ICL = ', results_wavsep[0], 'Flux ICL = ', results_wavsep[3] ) )
                print('SIZESEP | %12s%1.3e | %12s%1.3e | %12s%1.3e |  ' %('Flux gal = ', results_wavsep[7], 'Flux gal = ', results_wavsep[1], 'Flux gal = ', results_wavsep[4] ) )
                print('SIZESEP | %12s%1.3e | %12s%1.3e | %12s%1.3e | ' %('fICL = ', results_wavsep[8], 'fICL = ', results_wavsep[2], 'fICL = ', results_wavsep[5] ) )

            #-------------------------------------------------------------------
            # SBT
            for sbt in sbt_l:
                norm = header['NORM']
                make_results_sbt(oim, nfp, lvl_sep_big, sbt, norm, xs, ys, n_levels, plot_vignet = True)
                results_wavsep = measure_icl_quantities_sbt( oim, nfp, gamma = gamma, pixscale = pixscale, lvl_sep_big = lvl_sep_big, sbt = sbt, norm = norm, n_levels = n_levels, r_lsst = r_lsst, verbose = False  )
                print('SBT | %12s%7.1f |' %('mu = ', sbt ))
                print('SBT | %12s%1.3e |' %( 'Flux ICL = ', results_wavsep[0] ) )
                print('SBT | %12s%1.3e |  ' %( 'Flux gal = ', results_wavsep[1] ) )
                print('SBT | %12s%1.3e | ' %( 'fICL = ', results_wavsep[2] ) )

            #rdc, gal, iclbcg, res, rim = make_results_hardsepBCG( oim, nfp, lvl_sep_big, lvl_sep, rc, ricl, nf, xs, ys, n_levels )

            #print('HARDSEPBCG | LVL = %d ' %lvl_sep)
            #print('HARDSEPBCG | Flux ICL + BCG = %1.3e' %np.sum(iclbcg) )
            #print('HARDSEPBCG | Flux gal = %1.3e' %np.sum(gal) )
            #print('HARDSEPBCG | Fraction ICL + BCG = %1.3f' %( np.sum(iclbcg) / (np.sum(gal) + np.sum(iclbcg)) ) )
