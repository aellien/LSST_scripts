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
from glob import glob
from astropy.io import fits
from astropy.visualization import LinearStretch, LogStretch
from astropy.visualization import ZScaleInterval, MinMaxInterval
from astropy.visualization import ImageNormalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import SymLogNorm
from scipy.stats import sigmaclip

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def make_results( oim, path_wavelets, n_softhard_icl, n_hard_icl, rc, nf, xs, ys, n_levels ):

    # path, list & variables
    res = np.zeros( (xs, ys) )
    icl = np.zeros( (xs, ys) )
    gal = np.zeros( (xs, ys) )
    rim = np.zeros( (xs, ys) )

    rdc_array = np.zeros( (xs, ys, n_levels) )

    xc = xs / 2.
    yc = ys / 2.

    # Read atoms
    nf = nf[:-4]
    nf = ''.join( (nf, 'ol.*') )
    rimpath = os.path.join(path_wavelets, nf)

    for it in glob(rimpath):

        print(it)
        ol = d.read_objects_from_pickle( it )
        atom = np.zeros(oim.shape)

        for object in ol:


            x_min, y_min, x_max, y_max = object.bbox
            xco = x_min + ( x_max - x_min ) / 2
            yco = y_min + ( y_max - y_min ) / 2

            rdc_array[ x_min : x_max, y_min : y_max, object.level ] += object.image * gamma

            if object.level >= n_softhard_icl:

                if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) <= rc * 2**object.level:
                    icl[ x_min : x_max, y_min : y_max ] += object.image
                else:
                    gal[ x_min : x_max, y_min : y_max ] += object.image
            else:
                gal[ x_min : x_max, y_min : y_max ] += object.image * gamma

    # Datacube
    rdc = d.datacube( rdc_array, dctype = 'NONE', fheader = header )
    rim = np.sum( rdc_array, axis = 2 )
    res = oim - rim

    # write to fits
    hduo = fits.PrimaryHDU(res)
    hduo.writeto(os.path.join( path_wavelets, ''.join( ( nf, 'results.residuals.fits') )), overwrite = True )

    hduo = fits.PrimaryHDU(icl)
    hduo.writeto(os.path.join( path_wavelets, ''.join( ( nf, 'results.icl.fits') )), overwrite = True )

    hduo = fits.PrimaryHDU(gal)
    hduo.writeto(os.path.join( path_wavelets, ''.join( ( nf, 'results.gal.fits') )), overwrite = True )

    hduo = fits.PrimaryHDU(rim)
    hduo.writeto(os.path.join( path_wavelets, ''.join( ( nf, 'results.rim.fits') )), overwrite = True )

    return rdc, icl, gal, res, rim

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Plots

def plot_dawis_results( oim, oicl, ogal, rdc, icl, gal, res, rim, path_plots ):

    fig = plt.figure()
    ax = fig.subplots(3, 3, sharex = True, sharey = True, \
                                gridspec_kw = { 'hspace': 0.2, 'wspace': 0.2 })

    #--------------------------------------------------------------------
    # Original GAL+ICL
    axim = ax[0][0].imshow(oim,  norm = ImageNormalize( oim, \
                                            interval = ZScaleInterval(), \
                                            stretch = LinearStretch()), \
                                            origin = 'lower', \
                                            cmap = 'viridis' )
    divider = make_axes_locatable(ax[0][0])
    cax = divider.append_axes( 'right', size = "5%", pad = 0.05)
    caxor = fig.colorbar( axim, cax = cax, \
                                    orientation = 'vertical',
                                    pad = 0,
                                    shrink = 1.0 )
    caxor.set_label(r'ADU')
    ax[0][0].set_title('Original ICL + Galaxies')

    #--------------------------------------------------------------------
    # Original GAL
    axim = ax[0][1].imshow(ogal, norm = ImageNormalize( oim, \
                                      interval = ZScaleInterval(),
                                      stretch = LinearStretch()), \
                                      origin = 'lower', \
                                      cmap = 'viridis' )

    divider = make_axes_locatable(ax[0][1])
    cax = divider.append_axes( 'right', size = "5%", pad = 0.05)
    caxor = fig.colorbar( axim, cax = cax, \
                                    orientation = 'vertical',
                                    pad = 0,
                                    shrink = 1.0 )
    caxor.set_label(r'ADU')
    ax[0][1].set_title('Original Galaxies')

    #-------------------------------------------------------------------
    # Original ICL
    axim = ax[0][2].imshow(oicl, norm = ImageNormalize( oicl, \
                                      interval = MinMaxInterval(),
                                      stretch = LinearStretch()), \
                                      origin = 'lower', \
                                      cmap = 'viridis_r' )

    divider = make_axes_locatable(ax[0][2])
    cax = divider.append_axes( 'right', size = "5%", pad = 0.05)
    caxor = fig.colorbar( axim, cax = cax, \
                                    orientation = 'vertical',
                                    pad = 0,
                                    shrink = 1.0 )
    caxor.set_label(r'ADU')
    ax[0][2].set_title('Original ICL')

    #-------------------------------------------------------------------
    # Restored GAL+ICL
    axim = ax[1][0].imshow(rim, norm = ImageNormalize( oim, \
                                      interval = ZScaleInterval(),
                                      stretch = LinearStretch()), \
                                      origin = 'lower', \
                                      cmap = 'viridis' )

    divider = make_axes_locatable(ax[1][0])
    cax = divider.append_axes( 'right', size = "5%", pad = 0.05)
    caxor = fig.colorbar( axim, cax = cax, \
                                    orientation = 'vertical',
                                    pad = 0,
                                    shrink = 1.0 )
    caxor.set_label(r'ADU')
    ax[1][0].set_title('Restored ICL+Galaxies')

    #-------------------------------------------------------------------
    # Restored GAL
    axim = ax[1][1].imshow( gal, norm = ImageNormalize( oim, \
                                      interval = ZScaleInterval(),
                                      stretch = LinearStretch()), \
                                      origin = 'lower', \
                                      cmap = 'viridis' )

    divider = make_axes_locatable(ax[1][1])
    cax = divider.append_axes( 'right', size = "5%", pad = 0.05)
    caxor = fig.colorbar( axim, cax = cax, \
                                    orientation = 'vertical',
                                    pad = 0,
                                    shrink = 1.0 )
    caxor.set_label(r'ADU')
    ax[1][1].set_title('Restored Galaxies')

    #-------------------------------------------------------------------
    # Restored ICL
    axim = ax[1][2].imshow( icl, norm = ImageNormalize( oim - gal, \
                                      interval = ZScaleInterval(), \
                                      stretch = LinearStretch()), \
                                      origin = 'lower', \
                                      cmap = 'viridis' )

    divider = make_axes_locatable(ax[1][2])
    cax = divider.append_axes( 'right', size = "5%", pad = 0.05)
    caxor = fig.colorbar( axim, cax = cax, \
                                    orientation = 'vertical',
                                    pad = 0,
                                    shrink = 1.0 )
    caxor.set_label(r'ADU')
    ax[1][2].set_title('Restored ICL')

    #-------------------------------------------------------------------
    # Residuals GAL+ICL
    axim = ax[2][0].imshow(res,  norm = ImageNormalize( res, \
                                            interval = ZScaleInterval(), \
                                            stretch = LinearStretch()), \
                                            origin = 'lower', \
                                            cmap = 'viridis' )
    divider = make_axes_locatable(ax[2][0])
    cax = divider.append_axes( 'right', size = "5%", pad = 0.05)
    caxor = fig.colorbar( axim, cax = cax, \
                                    orientation = 'vertical',
                                    pad = 0,
                                    shrink = 1.0 )
    caxor.set_label(r'ADU')
    ax[2][0].set_title('Residuals ICL + Galaxies')

    #-------------------------------------------------------------------
    # Residuals GAL without ICL
    axim = ax[2][1].imshow(ogal - gal ,  norm = ImageNormalize( ogal - gal, \
                                            interval = ZScaleInterval(), \
                                            stretch = LinearStretch()), \
                                            origin = 'lower', \
                                            cmap = 'viridis' )
    divider = make_axes_locatable(ax[2][1])
    cax = divider.append_axes( 'right', size = "5%", pad = 0.05)
    caxor = fig.colorbar( axim, cax = cax, \
                                    orientation = 'vertical',
                                    pad = 0,
                                    shrink = 1.0 )
    caxor.set_label(r'ADU')
    ax[2][1].set_title('Residuals Galaxies without ICL')

    #-------------------------------------------------------------------
    # Residuals GAL leftover ICL
    axim = ax[2][2].imshow(oim - gal ,  norm = ImageNormalize( oim - gal, \
                                            interval = ZScaleInterval(), \
                                            stretch = LinearStretch()), \
                                            origin = 'lower', \
                                            cmap = 'viridis' )
    divider = make_axes_locatable(ax[2][2])
    cax = divider.append_axes( 'right', size = "5%", pad = 0.05)
    caxor = fig.colorbar( axim, cax = cax, \
                                    orientation = 'vertical',
                                    pad = 0,
                                    shrink = 1.0 )
    caxor.set_label(r'ADU')
    ax[2][2].set_title('Residuals Galaxies leftover ICL')

    '''
    axim = ax[2].imshow(nres, norm = SymLogNorm(1E1), vmin = -5E1, vmax = 5E1, \
                                            origin = 'lower', \
                                            cmap = 'seismic_r' )

    divider = make_axes_locatable(ax[2])
    cax = divider.append_axes( 'right', size = "5%", pad = 0.05)
    caxor = fig.colorbar( axim, cax = cax, \
                                    orientation = 'vertical',
                                    pad = 0,
                                    shrink = 1.0 )
    caxor.set_label(r'$\Delta$R[%]')
    ax[2].set_title('Residuals')
    '''

    plt.tight_layout()
    plt.show()
    #plt.savefig( os.path.join( path_plots, 'Euclid_dawis_results.pdf' ), format = 'pdf' )

if __name__ == '__main__':

    # Matplotlib params
    mpl.rcParams['xtick.major.size'] = 0
    mpl.rcParams['ytick.major.size'] = 0
    mpl.rcParams['xtick.minor.size'] = 0
    mpl.rcParams['ytick.minor.size'] = 0
    #mpl.rcParams['xtick.labelbottom'] = False
    #mpl.rcParams['ytick.labelleft'] = False

    # Paths, lists & variables
    path_data = '/home/ellien/LSST_ICL/simulations/out1'
    path_scripts = '/home/ellien/LSST_ICL/scripts'
    path_plots = '/home/ellien/LSST_ICL/plots/out1'
    path_wavelets = '/n03data/ellien/LSST_ICL/wavelets/out1/'
    gamma = 0.2
    n_levels = 11
    n_softhard_icl = 5
    n_hard_icl = 9
    rc = 100 # pixels, distance to center to be classified as ICL
    nfl = [ 'noise_00761_0000078_0.05_g_2Mpc.rebin.fits', \
            'noise_00761_0000120_0.05_g_2Mpc.rebin.fits', \
            'noise_00761_0000126_0.05_g_2Mpc.rebin.fits', \
            'noise_00761_0000174_0.05_g_2Mpc.rebin.fits', \
            'noise_00761_0000385_0.05_g_2Mpc.rebin.fits'  ]

    for nf in nfl:

        # Read files
        oimfile = os.path.join( path_data, nf )
        hdu = fits.open(oimfile)
        header = hdu[0].header
        oim = hdu[0].data

        xs, ys = oim.shape

        rdc, icl, gal, res, rim = make_results( oim, path_wavelets, n_softhard_icl, n_hard_icl, rc, nf, xs, ys, n_levels )
        #plot_dawis_results( oim, oicl, ogal, rdc, icl, gal, res, rim, path_plots )
