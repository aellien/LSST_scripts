#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# last modif : 06/2022
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import numpy as np
from astropy.io import fits
import os
import glob
from atom_props import *

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def measure_icl_quantities_sizesep( oim, nfp, gamma, lvl_sep_big, n_levels, n_iter = 1000, pdf = 'uniform', verbose = False ):

    # Paths, list & variables
    flux_icl_l   = []
    err_wr_icl_l = []
    flux_gal_l   = []
    err_wr_gal_l = []
    frac_icl_l   = []
    data_atom    = []

    # Read atoms and interscale trees
    ol, itl = read_image_atoms( nfp, verbose = True )

    # Compute mean size
    df_sizes = average_size_atom( ol, n_levels )

    if verbose:
        print(df_sizes)

    # Find wavelet scale of size derivative maximum
    lvl_sx = np.argmax( df_sizes['dsx_n'][2:] ) + 2
    lvl_sy = np.argmax( df_sizes['dsy_n'][2:] ) + 2

    lvl_sep = np.argmax(  ( df_sizes['dsx_n'][2:] + df_sizes['dsy_n'][2:] ) / 2  )

    # Get mean size of size derivative maximum
    sx = df_sizes['<sx>'][lvl_sep]
    low_sx = df_sizes['low<sx>'][lvl_sep]
    up_sx = df_sizes['up<sx>'][lvl_sep]

    if sx == 0.:
        sx = df_sizes['<sx>'][lvl_sep + 1]
        low_sx = df_sizes['low<sx>'][lvl_sep + 1]
        up_sx = df_sizes['up<sx>'][lvl_sep + 1]

    sy = df_sizes['<sy>'][lvl_sep]
    low_sy = df_sizes['low<sy>'][lvl_sep]
    up_sy = df_sizes['up<sy>'][lvl_sep]

    if sy == 0.:
        sy = df_sizes['<sx>'][lvl_sep +1]
        low_sy = df_sizes['low<sy>'][lvl_sep +1]
        up_sy = df_sizes['up<sy>'][lvl_sep +1]

    if verbose:
        print( 'x axis | level = %d  | <sx> = %d pix | low<sx> = %d | up<sx> = %d ' %(lvl_sx, sx, low_sx, up_sx))
        print( 'y axis | level = %d  | <sx> = %d pix | low<sy> = %d | up<sy> = %d ' %(lvl_sy, sy, low_sy, up_sy))

    # Store atom properties of interest in numpy array for speed up compuatations
    for j, o in enumerate(ol):
        x_min, y_min, x_max, y_max = o.bbox
        if o.level >= lvl_sep_big:
            data_atom.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
        else:
            data_atom.append( [ x_max - x_min, y_max - y_min, np.sum(o.image) * gamma, o.level ] )
    data_atom = np.array(data_atom)

    # Monte Carlo sampling
    for k in range( n_iter ):

        # Draw instance
        if pdf == 'normal':
            sig = 0.1
            ksx = np.random.normal( sx, sig * sx  )
            ksy = np.random.normal( sy, sig * sy  )
            print('Normal distribution, assuming std = %f x average size' %(sig))

        elif pdf == 'uniform':
            ksx = np.random.uniform( low_sx, up_sx  )
            ksy = np.random.uniform( low_sy, up_sy  )

        # Atom size separation
        x1 = np.where( data_atom[:,0] >= ksx )[0]
        x2 = np.where( data_atom[:,1] >= ksy )[0]
        xicl = np.unique( np.append( x1, x2 ) )
        mask = np.zeros( data_atom[:,0].shape, dtype = 'bool')
        mask[xicl] = True

        # ICL flux
        flux_icl = np.sum(data_atom[xicl][:,2])
        #err_wr_icl = np.sqrt( np.sum( afs[xicl][:, 5]**2 ))

        # Galaxy flux
        xgal = np.where(mask == False)[0]
        flux_gal = np.sum(data_atom[xgal][:,2])
        #err_wr_gal = np.sqrt( np.sum( afs[xgal][:, 5]**2 ))

        flux_icl_l.append( flux_icl )
        #err_wr_icl_l.append( err_wr_icl )
        flux_gal_l.append( flux_gal )
        #err_wr_gal_l.append( err_wr_gal )
        frac_icl_l.append( flux_icl / ( flux_icl + flux_gal) )

    return flux_icl_l, flux_gal_l, frac_icl_l, err_wr_icl_l, err_wr_gal_l

if __name__ == '__main__':

    # Paths, lists & variables
    path_data = '/home/ellien/LSST_ICL/simulations/out2/'
    path_scripts = '/home/ellien/LSST_ICL/scripts'
    path_wavelets = '/home/ellien/LSST_ICL/wavelets/out2/'

    dirl = ['HorizonAGN', 'Hydrangea', 'Magneticum' ] #, 'TNG-100']

    gamma = 0.5
    n_levels = 11
    lvl_sep_big = 6
    n_hard_icl = 6
    pixscale = 1.6 # ''/pixel
    physscale = 1 # kpc/''
    rc = 20 # pixels, distance to center to be classified as gal
    ricl = 5000 # pixels, distance to center to be classified as ICL

    '''
    for dir in dirl:

        print('\n\n%s --------------------------------------\n\n'%dir)

        image_dir = os.path.join( path_data, dir )
        image_files = glob.glob(image_dir+'/*.fits')
        image_files.sort()

        for nf in image_files:

            # Read files
            oimfile = os.path.join( path_data, nf )
            hdu = fits.open(oimfile)
            header = hdu[0].header
            oim = hdu[0].data
            norm = header['NORM']

            xs, ys = oim.shape

            split = nf.split('/')
            nf = split[-1]
            nfp = os.path.join( path_wavelets, dir, 'run1' )

            hdu = fits.open( os.path.join( nfp, nf[:-5] + '.results.icl.hardsep.fits' ))
            icl = hdu[0].data
            hicl = hdu[0].header

            hdu = fits.open( os.path.join( nfp, nf[:-5] + '.results.gal.hardsep.fits'  ))
            gal = hdu[0].data
            hgal = hdu[0].header

            hdu = fits.open( os.path.join( nfp, nf[:-5] + '.results.iclbcg.hardsepBCG.fits' ))
            iclbcg = hdu[0].data
            hicl = hdu[0].header

            hdu = fits.open( os.path.join( nfp, nf[:-5] + '.results.gal.hardsepBCG.fits'  ))
            galbcg = hdu[0].data
            hgal = hdu[0].header

            Ficl = np.sum(icl) * norm
            Fgal = np.sum(gal) * norm
            ficl = Ficl / (Ficl + Fgal)

            Ficlbcg = np.sum(iclbcg) * norm
            Fgalbcg = np.sum(galbcg) * norm
            ficlbcg = Ficlbcg / (Ficlbcg + Fgalbcg)

            print('\n%s'%nf)
            print('ICL:')
            print('Ficl = %e    Fgal = %e    ficl = %f'%(Ficl, Fgal, ficl))

            print('ICL + BCG:')
            print('Ficlbcg = %e    Fgal = %e    ficlbcg = %f'%(Ficlbcg, Fgalbcg, ficlbcg))
    '''
