import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
from scipy.optimize import curve_fit as cf
import os, random
#import fitting
import sys

#plt.style.use('mystyle') #path C:\Users\alexw\AppData\Local\Programs\Python\Python37\Lib\site-packages\matplotlib\mpl-data\stylelib

def findval(array,val):
    array = np.asarray(array)
    ind = np.abs(array - val).argmin()
    return ind

#imports the  FITTED spectra from the spectra folder
oldspecnames = next(os.walk('Spectra'))[2]
specnames = []
for i in [0,10,1000]:
    name = oldspecnames[i][0:20]+'-Prefitted.fits'
    specnames.append(name)

specstackinds = [0,1,2]

######-----stacking code----#######
fillval = np.nan
runlen = len(specstackinds)
highreslen = 50000
cutinds = []
normspecstore = np.empty([runlen,highreslen])
contspecstore = np.empty([runlen,highreslen])
normspecstore[:] = fillval
contspecstore[:] = fillval

wlenhighres = np.linspace(500, 4500, highreslen)
wlenmin = 10000
wlenmax = 0
stackstatus = []
specnumber = 0



for ind in specstackinds:
    specdata = fits.getdata('Fitted Spectra/'+specnames[ind],ext=1)#import fits image
    wlen = specdata.field(0)
    flux = specdata.field(1)
    contfit = specdata.field(6)
    specmetadata = fits.getdata('Fitted Spectra/'+specnames[ind],ext=2)#import fits image
    redshift = specmetadata[0][0]
    gcredshift = redshift+1

    cutstart = findval(wlen,(1030*(1+gcredshift)))
    cutend = findval(wlen,(1600*(1+gcredshift)))
    if cutend != 0 and cutstart != 0:
        cut = np.arange(cutstart,cutend)
        contfit = contfit[cut]
        flux = flux[cut]
        wlen = wlen[cut]

    wlenshift = wlen/(1+gcredshift) #+np.random.random(1)-0.5
    normspec = flux/contfit

    # plt.plot(wlenshift,flux)
    # plt.plot(wlenshift,contfit)
    #plt.plot(wlenshift,normspec)
    #plt.show()

    wlenintpol = interpolate.interp1d(wlenshift, normspec, 'linear', bounds_error=False, fill_value=fillval)
    contintpol = interpolate.interp1d(wlenshift, contfit, 'linear', bounds_error=False, fill_value=fillval)
    if wlenshift[0] < wlenmin:
        wlenmin = wlenshift[0]
    if wlenshift[-1] > wlenmax:
        wlenmax = wlenshift[-1]




    normspecstore[specnumber, 0:] = wlenintpol(wlenhighres)
    contspecstore[specnumber, 0:] = contintpol(wlenhighres)

    specnumber = specnumber + 1


#stacking data only if there are spec to stack

#cut extra zeropadding
start = (np.abs(wlenhighres - wlenmin)).argmin()# +10
end = (np.abs(wlenhighres - wlenmax)).argmin() #-10
wlenhighres = wlenhighres[start:end]
#cut downcols / remove empty rows
#normspecstore = np.delete(normspecstore, cutinds, axis = 0)
normspecstore = normspecstore[0:, start:end]
meanspec = np.nanmean(normspecstore, axis=0)
medspec = np.nanmedian(normspecstore, axis=0)

plt.plot(wlenhighres,normspecstore[0, 0:])
plt.plot(wlenhighres,normspecstore[1, 0:])
plt.plot(wlenhighres,normspecstore[2, 0:])
plt.show()
