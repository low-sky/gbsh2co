"""

Create a catalog of SCUBA2 clumps images for Band 5 follow up

"""

import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy import ndimage
from skimage.feature import blob_dog, blob_log, blob_doh
import scipy.ndimage as nd
import astropy.wcs as wcs
import sys


def box_mad(image, r):
    
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




hdu = fits.open(sys.argv[1])

data = hdu[0].getdata()
error = hdu[1].getdata()

# Chomp astrometry
w = wcs.WCS(fits.getheader(hdu[0].getheader())
            

blobs_doh = blob_doh(data, max_sigma=30, threshold=.01)
ra,dec = w.wcs_pix2sky(x,y,0)

