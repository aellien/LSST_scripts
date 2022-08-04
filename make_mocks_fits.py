from astropy.io import fits
import numpy as np
import os
import h5py
import glob
from skimage.transform import downscale_local_mean
from PIL import Image
import copy

def NormalizeData(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))

def make_img_map(img, mu_lim=32, pixel_size=0.2, rng=[22,36]):
    imgplot = img + np.random.normal(loc=0, scale=(10**(-0.4*mu_lim)*pixel_size*10.)/3., size=img.shape)
    imgflux = copy.copy(imgplot)
    imgplot = -2.5*np.log10(imgplot) + 2.5*np.log10((pixel_size)**2)
    imgplot[~np.isfinite(imgplot)] = rng[1]
    imgplot[imgplot < rng[0]] = rng[0]
    imgplot[imgplot > rng[1]] = rng[1]
    imgplot = NormalizeData(imgplot)
    imgplot = np.uint8(imgplot * 255)

    return imgplot, imgflux

if __name__ == '__main__':

    # Paths, lists & variables
    path_data = '/n03data/ellien/LSST_ICL/simulations/out2'
    path_scripts = '/home/ellien/LSST_ICL/scripts'
    path_output = '/n03data/ellien/LSST_ICL/simulations/out4'

    dirl = ['TNG-100']

    for dir in dirl:

        image_dir = os.path.join( path_data, dir )
        print(image_dir)
        image_files = glob.glob(image_dir+'/*001[57]*.hdf5')
        image_i = 0
        mu_lim = 30.3 # limiting SB, 3 sigma 10" X 10"

        make_png = False
        make_fits = True

        print(image_files)

        for image_file in image_files:

            print(image_file)
            with h5py.File(image_file, 'r') as f:
                image = f['image'][()]

            imgplot, imgflux = make_img_map(image, mu_lim=mu_lim, rng=[22,mu_lim+3])

            if make_png:
                imgname = image_files[0][:-5]+'.png'
                image_plot = Image.fromarray(np.uint8(downscale_local_mean(imgplot,(4,4))))
                image_plot.save(imgname)

            if make_fits:
                # produce fits image with scaling so that -2.5 log10(f) gives AB magnitudes
                fitsname = image_file[:-5]+'_mu%2.1f.fits' %(mu_lim)
                hdu = fits.PrimaryHDU(imgflux)
                hdul = fits.HDUList([hdu])
                hdul.writeto( os.path.join( image_dir, fitsname ), overwrite = True )

            print('done')
