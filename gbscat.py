import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy import ndimage
import skimage
from skimage.feature import blob_dog, blob_log, blob_doh, peak_local_max
import scipy.ndimage as nd
import astropy.wcs as wcs
import sys


def box_mad(image, r):
    """

    Create a catalog of SCUBA2 clumps images for Band 5 follow up.

    """    

    nd.filters.median_filter()


def MAD(a, c=0.6745, axis=None):
    """
    Adam Ginsburg's
    Median Absolute Deviation along given axis of an array:
    median(abs(a - median(a))) / c
    c = 0.6745 is the constant to convert from MAD to std; it is used by
    default
    """

    a = ma.masked_where(a!=a, a)
    if a.ndim == 1:
        d = ma.median(a)
        m = ma.median(ma.fabs(a - d) / c)
    else:
        d = ma.median(a, axis=axis)

    # I don't want the array to change so I have to copy it?
    if axis > 0:
        aswp = ma.swapaxes(a,0,axis)
    else:
        aswp = a
    m = ma.median(ma.fabs(aswp - d) / c, axis=0)
    return m



#file = sys.argv[1]
file = 'Aquila_450.fits'
hdu = fits.open(file)

data = np.squeeze(hdu[0].data)
var = np.squeeze(hdu[1].data)

# Chomp astrometry
w = wcs.WCS((hdu[0].header))

pc90  = np.percentile(var[np.isfinite(var)],90)
pc50 = np.percentile(var[np.isfinite(var)],50)
mask = (var < pc90)
circ = skimage.morphology.disk(15)
mask = ndimage.morphology.binary_erosion(mask,structure = circ)
mask = mask*(data>0)

data[~mask]=0.

pks = peak_local_max(data,min_distance=50,threshold_abs=3*np.sqrt(pc50))
y = pks[:,0]
x = pks[:,1]
#mxfilt = nd.maximum_filter(data,size=50,mode='constant')

#blobs_doh = blob_doh(data, max_sigma=50, threshold=np.sqrt(pc50)/10,
#                     min_sigma=5)
#x = blobs_doh[:,1]
#y = blobs_doh[:,0]

ra,dec,_ = w.wcs_pix2world(x,y,np.zeros(len(x)),0)

import matplotlib.pyplot as plt

plt.imshow(data,vmin=0,vmax=1e-1,cmap='copper')
plt.colorbar()
plt.scatter(x,y,marker='H',color='blue',facecolor='None')
plt.show()

