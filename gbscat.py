import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy import ndimage
import skimage
from skimage.feature import blob_dog, blob_log, blob_doh, peak_local_max
import scipy.ndimage as nd
import astropy.wcs as wcs
import sys
import matplotlib.pyplot as plt
import os
import subprocess

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

def h2cocat(filelist):
    for file in filelist:
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



filelist = np.loadtxt('filelist.txt',dtype='a')
ndfexec = '/home/ubuntu/star-hikianalia/bin/convert/ndf2fits'

t = Table(names=('RA','DEC','FLUX_450','XIMG','YIMG','FILENAME'),
          dtype=('f8','f8','f8','i4','i4','a40'))



for ndffile in filelist:
    base = os.path.basename(ndffile)
    parts = base.split('_')
    parts = parts[0:(np.where(np.char.strip(parts)=='450'))[0]]
    region = '_'.join(item for item in parts)
    file = region+'_450.fits'
    if not os.path.exists(file):
        callstring = ndfexec+' '+ndffile+' '+file
        subprocess.call(callstring,shell=True)

    #file = sys.argv[1]

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

    mask2 = data>(3*np.sqrt(var))
    circ = skimage.morphology.disk(1)
    mask2 = ndimage.morphology.binary_opening(mask2,structure = circ)

    mask = mask*mask2
    data[~mask]=0.

    pks = peak_local_max(data,min_distance=20,threshold_abs=3*np.sqrt(pc50))

    if len(pks)>0:

        y = pks[:,0]
        x = pks[:,1]

        fluxval = data[y,x]
        ra,dec,_ = w.wcs_pix2world(x,y,np.zeros(len(x)),0)

        for idx,_ in enumerate(x):
            t.add_row()
            t['XIMG'][-1]=x[idx]
            t['YIMG'][-1]=y[idx]
            t['RA'][-1]=ra[idx]
            t['DEC'][-1]=dec[idx]
            t['FLUX_450'][-1]=fluxval[idx]
            t['FILENAME'] = file

        realz = data[data>0]
        upperlim = np.percentile(realz,80)


        plt.imshow(data,vmin=0,vmax=upperlim,cmap='copper_r')
        plt.colorbar()
        plt.scatter(x,y,marker='H',color='blue',facecolor='cyan')
        #plt.scatter(x,y,marker='+',color='cyan',facecolor='None')
        plt.savefig(region+'.pdf')
        plt.close()
        plt.clf()
        t.write('s2_450.fits',overwrite=True,format='fits')

    #mxfilt = nd.maximum_filter(data,size=50,mode='constant')
    #blobs_doh = blob_doh(data, max_sigma=50, threshold=np.sqrt(pc50)/10,
    #                     min_sigma=5)
    #x = blobs_doh[:,1]
    #y = blobs_doh[:,0]
