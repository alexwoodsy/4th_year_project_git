import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
import os
#import fitting
import sys
#sys.path.append('C:/Users/jason/GIT/4th_year_project_git/Continuum Fitting')
sys.path.append('C:/Users/alexw/Documents/GitHub/4th_year_project_git/Continuum fitting')

def findval(array,val):
    array = np.asarray(array)
    ind = np.abs(array - val).argmin()
    return ind





#read in table with no signal to noise:
amdata = fits.getdata('Anti-Match/AM-nodups.fits',ext=1)#import fits image
amlen = len(amdata)
am_plate = amdata.field(3)
am_mjd = amdata.field(4)
am_fiberid = amdata.field(5)
am_z = amdata.field(6)

folderpath = 'E:/Sky Spectra Zip/Sky Spectra/'

stonall = np.zeros(amlen)
stonforest = np.zeros(amlen)

outcheck = 0
for ind in range(0,amlen):
    print(outcheck)
    plate = str(am_plate[ind]).zfill(4)
    mjd = str(am_mjd[ind]).zfill(5)
    fiberid = str(am_fiberid[ind]).zfill(4)
    redshift = am_z[ind]
    spec = 'spec-'+plate+'-'+mjd+'-'+fiberid+'.fits'
    ampath = folderpath+spec
    specdata = fits.getdata(ampath,ext=1)#import fits image

    wlen = 10**specdata.field(1)
    flux = specdata.field(0)
    ivar = specdata.field(2)
    ivar[ivar == 0] = 0.0000001

    lyalpha = 1215.67*(1+redshift)#calc lya using redshift #metadata
    lyalphaind = findval(wlen,lyalpha) #metadata



    std = (1/ivar)**0.5
    stonall[ind] = np.median(flux/std) #metadata
    stonforest[ind] = np.median(flux[0:lyalphaind]/std[0:lyalphaind])

    outcheck = outcheck + 1




hduold = fits.open('Anti-Match/AM-nodups.fits')
oldcols = hduold[1].columns


stonall_col = fits.Column(name='STON_ALL', array = stonall, format='F')
stonforest_col = fits.Column(name='STON_FOREST', array = stonforest, format='F')



metadata = fits.BinTableHDU.from_columns(oldcols+stonall_col+stonforest_col)

hduold.close()

primary = fits.PrimaryHDU()
hdul = fits.HDUList([primary, metadata])

hdul.writeto('Anti-Match/AM-sntest.fits',overwrite = True)
