#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# last modif : 06/2022
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import numpy as np
from astropy.io import fits
import os
import glob
from atom_props import *

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def measure_icl_quantities_sizesep( oim, nfp, gamma, lvl_sep_big, n_hard_icl, rc, ricl, nf, xs, ys, n_levels ):

    # path, list & variables
    res = np.zeros( (xs, ys) )

    rim = np.zeros( (xs, ys) )
    unclass = np.zeros( (xs, ys) )

    lres = []
    licl = []
    luicl = []
    lgal = []
    lrim = []
    lunclass = []

    rdc_array = np.zeros( (xs, ys, n_levels) )

    xc = xs / 2.
    yc = ys / 2.

    ol, itl = read_image_atoms( nfp, verbose = True )
    df_sizes = average_size_atom( ol, n_levels )
    print(df_sizes)

    lvl_sx = np.argmax( df_sizes['dsx_n'][2:] ) + 2
    lvl_sy = np.argmax( df_sizes['dsy_n'][2:] ) + 2

    sx = df_sizes['<sx>'][lvl_sx]
    low_sx = df_sizes['low<sx>'][lvl_sx]
    up_sx = df_sizes['up<sx>'][lvl_sx]

    if sx == 0.:
        sx = df_sizes['<sx>'][lvl_sx + 1]
        low_sx = df_sizes['low<sx>'][lvl_sx + 1]
        up_sx = df_sizes['up<sx>'][lvl_sx + 1]

    sy = df_sizes['<sy>'][lvl_sy]
    low_sy = df_sizes['low<sy>'][lvl_sy]
    up_sy = df_sizes['up<sy>'][lvl_sy]

    if sy == 0.:
        sy = df_sizes['<sx>'][lvl_sy +1]
        low_sy = df_sizes['low<sy>'][lvl_sy +1]
        up_sy = df_sizes['up<sy>'][lvl_sy +1]

    print(lvl_sx, sx, low_sx, up_sx)
    print(lvl_sy, sy, low_sy, up_sy)

    n_iter = 1000
    flux_icl_l = []
    flux_gal_l = []
    frac_icl_l = []
    for k in range( n_iter ):

        print(k)

        ksx = np.random.uniform( low_sx, up_sx  )
        ksy = np.random.uniform( low_sy, up_sy  )

        icl = np.zeros( (xs, ys) )
        gal = np.zeros( (xs, ys) )

        atom_gal = []
        coo_gal = []

        for j, object in enumerate(ol):

            x_min, y_min, x_max, y_max = object.bbox
            lvlo = object.level
            xco = itl[j].interscale_maximum.x_max
            yco = itl[j].interscale_maximum.y_max

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

        flux_icl_l.append( np.sum(icl) )
        flux_gal_l.append( np.sum(gal) )
        frac_icl_l.append( np.sum(icl) / ( np.sum(icl) + np.sum(gal) ) )
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

    return flux_icl_l, flux_gal_l, frac_icl_l

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
