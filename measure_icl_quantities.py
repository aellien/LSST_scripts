#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# last modif : 10/2022
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import numpy as np
from astropy.io import fits
import os
import glob
from atom_props import *
from make_results import *

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def create_circular_mask( h, w, center = None, radius = None ):

    if center is None: # use the middle of the image
        center = (int(w/2), int(h/2))
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return mask

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def measure_icl_quantities_sizesep( oim, nfp, gamma, size_sep, err_size, lvl_sep_big, n_levels, r_lsst, verbose = False ):

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

    mask = create_circular_mask( xs, ys, center = None, radius = r_lsst )
    flux_icl = np.sum( im_icl[mask] )
    flux_gal = np.sum( im_gal[mask] )
    frac_icl = flux_icl / ( flux_gal + flux_icl )

    # Up error computations
    atom_icl = []
    atom_gal = []
    im_icl = np.zeros((xs, ys))
    im_gal = np.zeros((xs, ys))

    for j, o in enumerate(ol):

        x_min, y_min, x_max, y_max = o.bbox
        sx = x_max - x_min
        sy = y_max - y_min

        if o.level >= lvl_sep_big:
            if (sx >= size_sep * ( 1 + err_size ) ) or (sy >= size_sep * ( 1 + err_size ) ):
                atom_icl.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
                im_icl[ x_min : x_max, y_min : y_max ] += o.image
            else:
                atom_gal.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
                im_gal[ x_min : x_max, y_min : y_max ] += o.image
        else:
            if (sx >= size_sep * ( 1 + err_size ) ) or (sy >= size_sep * ( 1 + err_size ) ):
                atom_icl.append( [ x_max - x_min, y_max - y_min, np.sum(o.image) * gamma, o.level ] )
                im_icl[ x_min : x_max, y_min : y_max ] += o.image * gamma
            else:
                atom_gal.append( [ x_max - x_min, y_max - y_min, np.sum(o.image) * gamma, o.level ] )
                im_gal[ x_min : x_max, y_min : y_max ] += o.image * gamma

    atom_gal = np.array(atom_gal)
    atom_icl = np.array(atom_icl)
    flux_up_icl = np.sum( im_icl[mask] )
    flux_up_gal = np.sum( im_gal[mask] )
    frac_up_icl = flux_up_icl / ( flux_up_gal + flux_up_icl )

    # Low error computations
    atom_icl = []
    atom_gal = []
    im_icl = np.zeros((xs, ys))
    im_gal = np.zeros((xs, ys))

    for j, o in enumerate(ol):

        x_min, y_min, x_max, y_max = o.bbox
        sx = x_max - x_min
        sy = y_max - y_min

        if o.level >= lvl_sep_big:
            if (sx >= size_sep * ( 1 + err_size ) ) or (sy >= size_sep * ( 1 - err_size ) ):
                atom_icl.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
                im_icl[ x_min : x_max, y_min : y_max ] += o.image
            else:
                atom_gal.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
                im_gal[ x_min : x_max, y_min : y_max ] += o.image
        else:
            if (sx >= size_sep * ( 1 + err_size ) ) or (sy >= size_sep * ( 1 - err_size ) ):
                atom_icl.append( [ x_max - x_min, y_max - y_min, np.sum(o.image) * gamma, o.level ] )
                im_icl[ x_min : x_max, y_min : y_max ] += o.image * gamma
            else:
                atom_gal.append( [ x_max - x_min, y_max - y_min, np.sum(o.image) * gamma, o.level ] )
                im_gal[ x_min : x_max, y_min : y_max ] += o.image * gamma

    atom_gal = np.array(atom_gal)
    atom_icl = np.array(atom_icl)
    flux_low_icl = np.sum( im_icl[mask] )
    flux_low_gal = np.sum( im_gal[mask] )
    frac_low_icl = flux_low_icl / ( flux_low_gal + flux_low_icl )

    return flux_icl, flux_gal, frac_icl, flux_up_icl, flux_up_gal, frac_up_icl, flux_low_icl, flux_low_gal, frac_low_icl

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def measure_icl_quantities_sizesep_bckup_test( oim, nfp, gamma, lvl_sep_big, n_levels, n_iter = 1000, pdf = 'uniform', verbose = False ):

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

    lvl_sep = np.argmax(  ( df_sizes['dsx_n'][2:] + df_sizes['dsy_n'][2:] ) / 2  ) + 2

    # Get mean size of size derivative maximum
    sx = df_sizes['<sx>'][lvl_sep]
    low_sx = df_sizes['low<sx>'][lvl_sep]
    up_sx = df_sizes['up<sx>'][lvl_sep]

    sy = df_sizes['<sy>'][lvl_sep]
    low_sy = df_sizes['low<sy>'][lvl_sep]
    up_sy = df_sizes['up<sy>'][lvl_sep]

    if verbose:
        print( 'x axis | level = %d  | <sx> = %d pix | low<sx> = %d | up<sx> = %d ' %(lvl_sep, sx, low_sx, up_sx))
        print( 'y axis | level = %d  | <sx> = %d pix | low<sy> = %d | up<sy> = %d ' %(lvl_sep, sy, low_sy, up_sy))

    # Store atom properties of interest in numpy array for speed up computations
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def measure_icl_quantities_wavsep( oim, nfp, gamma, lvl_sep_big, lvl_sep, n_levels, r_lsst, verbose = False ):

    # Paths, list & variables
    atom_icl = []
    atom_gal = []
    xs, ys = oim.shape
    im_icl = np.zeros((xs, ys))
    im_gal = np.zeros((xs, ys))

    # Read atoms and interscale trees
    ol, itl = read_image_atoms( nfp, verbose = True )

    # Atom wavelet scale separation
    for j, o in enumerate(ol):
        x_min, y_min, x_max, y_max = o.bbox
        if o.level >= lvl_sep_big:
            if o.level >= lvl_sep:
                atom_icl.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
                im_icl[ x_min : x_max, y_min : y_max ] += o.image
            else:
                atom_gal.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
                im_gal[ x_min : x_max, y_min : y_max ] += o.image
        else:
            if o.level >= lvl_sep:
                atom_icl.append( [ x_max - x_min, y_max - y_min, np.sum(o.image) * gamma, o.level ] )
                im_icl[ x_min : x_max, y_min : y_max ] += o.image * gamma
            else:
                atom_gal.append( [ x_max - x_min, y_max - y_min, np.sum(o.image) * gamma, o.level ] )
                im_gal[ x_min : x_max, y_min : y_max ] += o.image * gamma

    atom_gal = np.array(atom_gal)
    atom_icl = np.array(atom_icl)

    mask = create_circular_mask( xs, ys, center = None, radius = r_lsst )
    flux_icl = np.sum( im_icl[mask] )
    flux_gal = np.sum( im_gal[mask] )
    frac_icl = flux_icl / ( flux_gal + flux_icl )

    # Up error computations
    atom_icl = []
    atom_gal = []
    im_icl = np.zeros((xs, ys))
    im_gal = np.zeros((xs, ys))

    for j, o in enumerate(ol):
        x_min, y_min, x_max, y_max = o.bbox
        if o.level >= lvl_sep_big:
            if o.level >= lvl_sep - 1:
                atom_icl.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
                im_icl[ x_min : x_max, y_min : y_max ] += o.image
            else:
                atom_gal.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
                im_icl[ x_min : x_max, y_min : y_max ] += o.image
        else:
            if o.level >= lvl_sep - 1:
                atom_icl.append( [ x_max - x_min, y_max - y_min, np.sum(o.image) * gamma, o.level ] )
                im_icl[ x_min : x_max, y_min : y_max ] += o.image * gamma
            else:
                atom_gal.append( [ x_max - x_min, y_max - y_min, np.sum(o.image) * gamma, o.level ] )
                im_gal[ x_min : x_max, y_min : y_max ] += o.image * gamma

    atom_gal = np.array(atom_gal)
    atom_icl = np.array(atom_icl)
    flux_up_icl = np.sum( im_icl[mask] )
    flux_up_gal = np.sum( im_gal[mask] )
    frac_up_icl = flux_up_icl / ( flux_up_gal + flux_up_icl )

    # Low error computations
    atom_icl = []
    atom_gal = []
    im_icl = np.zeros((xs, ys))
    im_gal = np.zeros((xs, ys))

    for j, o in enumerate(ol):
        x_min, y_min, x_max, y_max = o.bbox
        if o.level >= lvl_sep_big:
            if o.level >= lvl_sep + 1:
                atom_icl.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
                im_icl[ x_min : x_max, y_min : y_max ] += o.image
            else:
                atom_gal.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
                im_gal[ x_min : x_max, y_min : y_max ] += o.image
        else:
            if o.level >= lvl_sep + 1:
                atom_icl.append( [ x_max - x_min, y_max - y_min, np.sum(o.image) * gamma, o.level ] )
                im_icl[ x_min : x_max, y_min : y_max ] += o.image * gamma
            else:
                atom_gal.append( [ x_max - x_min, y_max - y_min, np.sum(o.image) * gamma, o.level ] )
                im_gal[ x_min : x_max, y_min : y_max ] += o.image * gamma

    atom_gal = np.array(atom_gal)
    atom_icl = np.array(atom_icl)
    flux_low_icl = np.sum( im_icl[mask] )
    flux_low_gal = np.sum( im_gal[mask] )
    frac_low_icl = flux_low_icl / ( flux_low_gal + flux_low_icl )

    return flux_icl, flux_gal, frac_icl, flux_up_icl, flux_up_gal, frac_up_icl, flux_low_icl, flux_low_gal, frac_low_icl

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def measure_icl_quantities_sbt( oim, nfp, gamma, pixscale, lvl_sep_big, sbt, norm, n_levels, r_lsst, verbose = False  ):

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

    mask = create_circular_mask( xs, ys, center = None, radius = r_lsst )
    im_tot[~mask] = 0.
    im_tot *= norm # renormalize for SB
    im_tot[im_tot < 10E-30] = 10E-30 # get rid of nul pixels
    im_tot_sb = - 2.5 * np.log10(im_tot / pixscale**2 )

    pos_icl = np.where(im_tot_sb > sbt)
    flux_icl = np.sum( im_tot[pos_icl] )
    pos_gal = np.where(im_tot_sb <= sbt)
    flux_gal = np.sum( im_tot[pos_gal] )
    frac_icl = flux_icl / ( flux_gal + flux_icl )

    return flux_icl, flux_gal, frac_icl

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def measure_icl_quantities_spawavsep( oim, nfp, gamma, lvl_sep_big, cat_gal, rc, n_sig_gal, lvl_sep, xs, ys, n_levels, r_lsst, verbose = False ):

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

    for j, object in enumerate(ol):

        x_min, y_min, x_max, y_max = object.bbox
        lvlo = object.level
        xco = itl[j].interscale_maximum.x_max
        yco = itl[j].interscale_maximum.y_max

        if object.level >= lvl_sep_big:

            if object.level < lvl_sep:

                flag_gal = False
                for pos in cat_gal:

                    if np.sqrt( ( xco - pos[1] )**2 + ( yco - pos[0] )**2 ) <= rc:
                        flag_gal = True
                        break

                if flag_gal == True:
                    if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) <= rc: #BCG
                        icl[ x_min : x_max, y_min : y_max ] += object.image
                    else:
                        gal[ x_min : x_max, y_min : y_max ] += object.image

                if flag_gal == False:
                    icl[ x_min : x_max, y_min : y_max ] += object.image

            else:
                icl[ x_min : x_max, y_min : y_max ] += object.image

        else:

            if object.level < lvl_sep:

                flag_gal = False
                for pos in cat_gal:

                    if np.sqrt( ( xco - pos[1] )**2 + ( yco - pos[0] )**2 ) <= rc:
                        flag_gal = True
                        break

                if flag_gal == True:
                    if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) <= rc:
                        icl[ x_min : x_max, y_min : y_max ] += object.image * gamma
                    else:
                        gal[ x_min : x_max, y_min : y_max ] += object.image * gamma

                if flag_gal == False:
                    icl[ x_min : x_max, y_min : y_max ] += object.image * gamma

            else:
                icl[ x_min : x_max, y_min : y_max ] += object.image * gamma

    mask = create_circular_mask( xs, ys, center = None, radius = r_lsst )
    flux_icl = np.sum( icl[mask] )
    flux_gal = np.sum( gal[mask] )
    frac_icl = flux_icl / ( flux_gal + flux_icl )

    return flux_icl, flux_gal, frac_icl

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def measure_icl_quantities_bcgwavsep( oim, nfp, gamma, lvl_sep_big, rc, lvl_sep, xs, ys, n_levels, r_lsst, verbose = False ):

    # Paths, list & variables
    atom_icl = []
    atom_gal = []
    xs, ys = oim.shape
    im_icl = np.zeros((xs, ys))
    im_gal = np.zeros((xs, ys))
    xc = xs / 2.
    yc = ys / 2.
    # Read atoms and interscale trees
    ol, itl = read_image_atoms( nfp, verbose = True )
    for j, object in enumerate(ol):

        x_min, y_min, x_max, y_max = object.bbox
        lvlo = object.level
        xco = itl[j].interscale_maximum.x_max
        yco = itl[j].interscale_maximum.y_max

        if object.level >= lvl_sep_big:

            if object.level >= lvl_sep:
                im_icl[ x_min : x_max, y_min : y_max ] += object.image
            else:
                if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) <= rc:
                    im_icl[ x_min : x_max, y_min : y_max ] += object.image
                else:
                    im_gal[ x_min : x_max, y_min : y_max ] += object.image

        else:

            if object.level >= lvl_sep:
                im_icl[ x_min : x_max, y_min : y_max ] += object.image * gamma
            else:
                if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) <= rc:
                    im_icl[ x_min : x_max, y_min : y_max ] += object.image * gamma
                else:
                    im_gal[ x_min : x_max, y_min : y_max ] += object.image * gamma

    mask = create_circular_mask( xs, ys, center = None, radius = r_lsst )
    flux_icl = np.sum( im_icl[mask] )
    flux_gal = np.sum( im_gal[mask] )
    frac_icl = flux_icl / ( flux_gal + flux_icl )

    return flux_icl, flux_gal, frac_icl

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def measure_icl_quantities_bcgsizesep( oim, nfp, gamma, lvl_sep_big, rc, size_sep, xs, ys, n_levels, r_lsst, verbose = False ):

    # Paths, list & variables
    atom_icl = []
    atom_gal = []
    xs, ys = oim.shape
    im_icl = np.zeros((xs, ys))
    im_gal = np.zeros((xs, ys))
    xc = xs / 2.
    yc = ys / 2.
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
                if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) <= rc: #BCG
                    im_icl[ x_min : x_max, y_min : y_max ] += o.image
                else:
                    atom_gal.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
                    im_gal[ x_min : x_max, y_min : y_max ] += o.image
        else:
            if (sx >= size_sep) or (sy >= size_sep):
                atom_icl.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
                im_icl[ x_min : x_max, y_min : y_max ] += o.image * gamma
            else:
                if np.sqrt( ( xco - xc )**2 + ( yco - yc )**2 ) <= rc: #BCG
                    im_icl[ x_min : x_max, y_min : y_max ] += o.image * gamma
                else:
                    atom_gal.append( [ x_max - x_min, y_max - y_min, np.sum(o.image), o.level ] )
                    im_gal[ x_min : x_max, y_min : y_max ] += o.image * gamma

    mask = create_circular_mask( xs, ys, center = None, radius = r_lsst )
    flux_icl = np.sum( im_icl[mask] )
    flux_gal = np.sum( im_gal[mask] )
    frac_icl = flux_icl / ( flux_gal + flux_icl )

    return flux_icl, flux_gal, frac_icl


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
