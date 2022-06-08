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
from scipy.stats import sigmaclip, median_absolute_deviation
from skimage.measure import label, regionprops
from sklearn.utils import resample

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def bootstrap_error( data, n_iter, alpha = 0.95  ):

    n_sample = int( 0.75 * np.size( data ) )
    statfrac = np.zeros( (n_iter) )

    for i in range( 0, n_iter ):

        newsample = resample( data, n_samples = n_sample )
        statfrac[i] = np.mean( newsample )

    statfrac = np.sort(statfrac)
    plow = ( (1.- alpha) / 2. ) * 100.
    low = np.percentile( statfrac, plow )

    pup = ( alpha + (1. - alpha) / 2.) * 100
    up = np.percentile( statfrac, pup )

    return low, up

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def read_image_atoms( nfp, verbose = False ):

    # Object lists
    opath = nfp + '*ol*pkl'
    opathl = glob.glob(opath)
    opathl.sort()

    # Interscale tree lists
    itpath = nfp + '*itl*pkl'
    itpathl = glob.glob(itpath)
    itpathl.sort()

    tol = []
    titl = []

    if verbose:
        print('Reading %s.'%(opath))
        print('Reading %s.'%(itpath))

    for i, ( op, itlp ) in enumerate( zip( opathl, itpathl )):

        if verbose :
            print('Iteration %d' %(i), end ='\r')

        ol = d.read_objects_from_pickle( op )
        itl = d.read_interscale_trees_from_pickle( itlp )

        for j, o in enumerate(ol):

            tol.append(o)
            titl.append(itl[j])

    return tol, titl

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def average_size_atom( ol, n_levels ):

    sx = list( np.zeros( ( n_levels, 1 ) ))
    sy = list( np.zeros( ( n_levels, 1 ) ))

    for o in ol:
        x_min, y_min, x_max, y_max = o.bbox
        xs = x_max - x_min
        ys = y_max - y_min
        sx[o.level] = np.append( sx[o.level], xs )
        sy[o.level] = np.append( sy[o.level], ys )

    data = []
    for i in range(n_levels):
        lowx, upx = bootstrap_error( sx[i], 1000, alpha = 0.95  )
        lowy, upy = bootstrap_error( sy[i], 1000, alpha = 0.95  )
        data.append( [ i, 2**i, np.mean(sx[i]), np.mean(sy[i]), lowx, upx, lowy, upy, np.std(sx[i]), \
                        np.std(sy[i]), np.median(sx[i]), np.median(sy[i]), \
                        median_absolute_deviation(sx[i]), median_absolute_deviation(sy[i]), np.mean(sx[i] / 2**i), \
                        np.mean(sy[i] / 2**i), np.std(sx[i] / 2**i), np.std(sy[i] / 2**i) ] )

    for i in range(n_levels - 1):
        data[i].append( (data[i+1][2] - data[i][2]) / 2**i )
        data[i].append( (data[i+1][3] - data[i][3]) / 2**i )

    df = pd.DataFrame( data, columns = [ 'z', 'sz', '<sx>', '<sy>', 'low<sx>', 'up<sx>', 'low<sy>', 'up<sy>', 'std(sx)', 'std(sy)', \
                                        'med(sx)', 'med(sy)', 'mad(sx)', 'mad(sy)', '<sx_n>', \
                                        '<sy_n>', 'std(sx_n)', 'std(sy_n)', 'dsx_n', 'dsy_n'] )

    return df

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def find_minimum_size_z( sizes, start, end ):

    sl = sizes[start:end]
    min =  sizes[start]
    lvl = start
    flag = False

    for i, s in enumerate(sl):

        if s <= min and s != 0.:
            min = s
            lvl = i + start
            flag = True

        elif flag == True:
            break

        else:
            continue

    if lvl == start:
        print('WARING: z_minimum = z_start !')

    return lvl, min

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def derivative( sizes, start, end ):

    ds = []

    for i in range(start, end - 1):

        if ( sizes[i] == 0. ) and ( sizes[ i + 1] == 0. ):
            ds.append(0.)

        elif ( sizes[i] == 0. ) and ( sizes[ i + 1] != 0. ):
            if ( i >= 1 ) and ( sizes[i - 1] != 0. ):
                ds.append( ( sizes[i + 1] - sizes[i - 1] ) / 2**i )
            else:
                ds.append(0.)

        elif ( sizes[i] != 0. ) and ( sizes[ i + 1] == 0. ):
            if ( i < end - 2 ) and ( sizes[i + 2] != 0. ):
                ds.append( ( sizes[i + 2] - sizes[i] ) / 2**i )
            else:
                ds.append(0.)

        else:
            ds.append( ( sizes[i+1] - sizes[i] ) / 2**i )

    ds = np.array( ds )

    return ds


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def average_lum_atom( ol, n_levels ):

    sx = list( np.zeros( ( n_levels, 1 ) ))
    sy = list( np.zeros( ( n_levels, 1 ) ))

    for o in ol:
        x_min, y_min, x_max, y_max = o.bbox
        xs = x_max - x_min
        ys = y_max - y_min
        sx[o.level] = np.append( sx[o.level], xs )
        sy[o.level] = np.append( sy[o.level], ys )

    data = []
    for i in range(n_levels):
        data.append( [ i, 2**i, np.mean(sx[i]), np.mean(sy[i]), np.std(sx[i]), np.std(sy[i]), np.median(sx[i]), np.median(sy[i]), median_absolute_deviation(sx[i]), median_absolute_deviation([sy[i]]) ] )

    for i in range(n_levels - 1):
        data[i].append( (data[i+1][2] - data[i][2]) / 2**i )
        data[i].append( (data[i+1][3] - data[i][3]) / 2**i )

    df = pd.DataFrame( data, columns = [ 'z', 'sz', '<sx>', '<sy>', 'std(sx)', 'std(sy)', 'med(sx)', 'med(sy)', 'mad(sx)', 'mad(sy)', 'dsx_n', 'dsy_n'] )

    return df

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if __name__ == '__main__':

    # Paths, lists & variables
    path_data = '/n03data/ellien/LSST_ICL/simulations/out2/'
    path_scripts = '/home/ellien/LSST_ICL/LSST_scripts'
    path_wavelets = '/n03data/ellien/LSST_ICL/wavelets/out2/'
    path_plots = '/home/ellien/LSST_ICL/plots'

    dirl = ['HorizonAGN', 'Hydrangea', 'Magneticum', 'TNG-100']

    gamma = 0.5
    n_levels = 11
    lvl_sep_big = 5
    n_hard_icl = 5
    pixscale = 1.6 # ''/pixel
    physscale = 1 # kpc/''

    dfl = []
    data = []

    for dir in dirl:

        image_dir = os.path.join( path_data, dir )
        image_files = glob.glob(image_dir+'/*norm.fits')

        for nf in image_files:

            print('\n%s'%nf)

            # Read atoms
            nf = nf[:-4]
            split = nf.split('/')
            nf = split[-1]
            nfp = os.path.join( path_wavelets, dir, 'run1', nf )
            opath = nfp + '*ol*pkl'
            opathl = glob.glob(opath)
            opathl.sort()
            itpath = nfp + '*itl*pkl'
            itpathl = glob.glob(itpath)
            itpathl.sort()
            tol = []
            titl = []

            for i, ( op, itlp ) in enumerate( zip( opathl, itpathl )):

                print('read iteration %d' %(i), end ='\r')

                ol = d.read_objects_from_pickle( op )
                itl = d.read_interscale_trees_from_pickle( itlp )

                for j, o in enumerate(ol):

                    x_min, y_min, x_max, y_max = o.bbox
                    xs = x_max - x_min
                    ys = y_max - y_min
                    xmax = np.max(o.image)
                    xmean = np.mean(o.image)

                    tol.append(o)
                    titl.append(itl[j])
                    data.append( [ o.level, xs, ys, xmax, xmean ] )

            print('\n%d sources in total.' %len(tol))

    df = pd.DataFrame( data, columns = [ 'z', 'sx', 'sy', 'xmax', 'xmean' ] )
    mean_sx = []
    mean_sy = []
    norm_mean_sx = []
    norm_mean_sy = []
    mean_xmax = []
    mean_xmean = []
    std_sx = []
    std_sy = []
    norm_std_sx = []
    norm_std_sy = []
    std_xmax = []
    std_xmean = []
    lvls = []
    for i in range(n_levels):
        lvls.append(i)
        mean_sx.append( np.mean( df[ df['z'] == i ]['sx'] ) )
        mean_sy.append( np.mean( df[ df['z'] == i ]['sy'] ) )
        mean_xmax.append( np.median( df[ df['z'] == i ]['xmax'] ) )
        mean_xmean.append( np.median( df[ df['z'] == i ]['xmean'] ) )
        norm_mean_sx.append( np.mean( df[ df['z'] == i ]['sx'] / 2**df[ df['z'] == i ]['z'] ) )
        norm_mean_sy.append( np.mean( df[ df['z'] == i ]['sy'] / 2**df[ df['z'] == i ]['z'] ) )
        norm_std_sx.append( np.std( df[ df['z'] == i ]['sx'] / 2**df[ df['z'] == i ]['z'] ) )
        norm_std_sy.append( np.std( df[ df['z'] == i ]['sy'] / 2**df[ df['z'] == i ]['z'] ) )
        std_sx.append( np.std( df[ df['z'] == i ]['sx'] ) )
        std_sy.append( np.std( df[ df['z'] == i ]['sy'] ) )
        std_xmax.append( median_absolute_deviation( df[ df['z'] == i ]['xmax'] ) )
        std_xmean.append( median_absolute_deviation( df[ df['z'] == i ]['xmean'] ) )

    fig, ax = plt.subplots( 2, 2 )
    lvls = np.array(lvls)
    ax[0][0].errorbar( lvls, mean_sx, yerr = std_sx, color = 'blue', alpha = 0.5 )
    ax[0][0].errorbar( lvls, mean_sy, yerr = std_sy, color = 'red', alpha = 0.5 )

    ax[0][1].errorbar( lvls, mean_xmax, yerr = std_xmax, color = 'green', alpha = 0.5 )
    ax[0][1].set_yscale('log')

    ax[1][1].errorbar( lvls, norm_mean_sx, yerr = norm_std_sx, color = 'black', alpha = 0.5 )
    ax[1][1].errorbar( lvls, norm_mean_sy, yerr = norm_std_sy, color = 'red', alpha = 0.5 )

    ax[1][0].errorbar( lvls, mean_xmean, yerr = std_xmean, color = 'pink', alpha = 0.5 )
    ax[1][0].set_yscale('log')

    fig.savefig(os.path.join(path_plots, 'average_size_vs_z.pdf'), format = 'pdf')

    plt.figure()
    for i in range(3, n_levels):
        print(len(df[ df['z'] == i ]))
        df[ df['z'] == i ]['sx'].hist(bins = 100)
    plt.yscale('log')
    plt.xscale('log')
    plt.savefig(os.path.join(path_plots, 'average_size_hist_per_z.pdf'), format = 'pdf')

    plt.figure()
    plt.plot( df[df['z']>=4]['sx'], df[df['z']>=4]['xmax'], 'bo', markersize = 1, alpha = 0.05)
    plt.yscale('log')
    plt.xscale('log')
    plt.savefig(os.path.join(path_plots, 'average_size_vs_xmax.pdf'), format = 'pdf')

    plt.figure()
    plt.plot( df[df['z']>=4]['sx'], df[df['z']>=4]['xmean'], 'bo', markersize = 1, alpha = 0.05)
    plt.yscale('log')
    plt.xscale('log')
    plt.savefig(os.path.join(path_plots, 'average_size_vs_xmean.pdf'), format = 'pdf')

    plt.close()
