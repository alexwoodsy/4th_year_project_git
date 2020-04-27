import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate
import os

specnames = next(os.walk('Spectra'))[2] #dir is your directory path as string
spectot = len(specnames)
#add indexing for epctra in file to allow loop over all
#specind=1
specsample = np.array([0])

for specind in specsample:
    specdirectory = 'Spectra/'+specnames[specind]
    #print(specdirectory)
    data = fits.getdata(specdirectory,ext=1)#import fits image
    speclen = len(data)
    flux = np.zeros(speclen)
    wlen = np.zeros(speclen)

    for i in range(0,speclen):
     flux[i] = data[i][0]
     wlen[i] = 10**(data[i][1])


#meta data extraction to get z:
    fitdata = fits.getdata(specdirectory,ext=2)#import fits image
    metasize = len(fitdata[0])
    #print(metasize)
    if metasize == 126:
         redshift = fitdata[0][63]
    else:
         redshift = fitdata[0][38]


    lyalphacalc = 1215.67*(1+redshift)#calc lya using redshift
    print('lyalpha = ',lyalphacalc)
    lyalphaind = (np.abs(wlen - lyalphacalc)).argmin()#finds index of nearest point in data

#fitting:
    intervals = 15
    window = int(speclen/intervals)
    windowwlen = wlen[window]-wlen[0]

    step = 0
    i = 0
    intervalwlen = np.zeros(intervals+2)
    winpeak = np.zeros(intervals+2)

    def findmax(array):
        array = np.asarray(array)
        ind = (np.abs(array - np.max(array))).argmin()
        return ind

    def findmed(array):
        array = np.asarray(array)
        ind = (np.abs(array - np.median(array))).argmin()
        return ind

    while (step+window) <= speclen:
        windata = flux[step:(step+window)]
        winpeakmed = step + findmed(windata)
        winpeakmax = step + findmax(windata)
        if wlen[winpeakmax] < wlen[lyalphaind]:
            winpeakind = winpeakmax
        elif wlen[winpeakmax] > wlen[lyalphaind] and wlen[winpeakmed] < wlen[lyalphaind]:
            winpeakind = winpeakmax
        else:
            winpeakind = winpeakmed

        winpeak[i+1] = flux[winpeakind]

#stops slection of interval near lyalpha peak
        if np.abs(wlen[winpeakind] - wlen[lyalphaind]) > window:
            intervalwlen[i+1] = wlen[winpeakind]
        else:
            intervalwlen[i+1] = winpeak[i+1] = 0

        step = step + window
        i = i + 1

#remove zero values made by filtering procedure
    winpeak = winpeak[winpeak != 0]
    intervalwlen = intervalwlen[intervalwlen != 0]
#pad interval with start/end value to allign correctly
    winpeak[0],winpeak[-1] = flux[0],flux[-1]
    intervalwlen[0],intervalwlen[-1] = wlen[0],wlen[-1]

    intpol = interpolate.interp1d(intervalwlen, winpeak, kind=1)
    contfit = intpol(wlen)
    normspec = flux-contfit


#plotting:
    wlim = speclen
    plt.plot(wlen[0:wlim],flux[0:wlim],label=specnames[specind])
    plt.plot(intervalwlen[0:wlim],winpeak[0:wlim],'*',label='intervals')
    plt.plot(wlen[0:wlim],contfit[0:wlim],'--',label='interpolation')
    plt.plot(wlen[lyalphaind],flux[lyalphaind],'.',label='lyalpha')
    plt.xlabel('wavelength (Angstroms)')
    plt.ylabel('Flux')
plt.legend()

plt.figure()
plt.plot(wlen[0:wlim],normspec[0:wlim],label='contrem')
plt.xlabel('wavelength (Angstroms)')
plt.ylabel('Flux')
plt.legend()
plt.show()
