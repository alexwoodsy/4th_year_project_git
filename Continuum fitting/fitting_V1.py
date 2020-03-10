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
    print(specdirectory)
    data = fits.getdata(specdirectory,ext=1)#import fits image
    speclen = len(data)
    flux = np.zeros(speclen)
    wlen = np.zeros(speclen)
    model = np.zeros(speclen)
    ivar = np.zeros(speclen)

    for i in range(0,speclen):
     flux[i] = data[i][0]
     ivar[i] = data[i][2]
     model[i] = data[i][7]
     wlen[i] = 10**(data[i][1])

    #ston calculation
    std = (ivar)**-0.5
    ston = np.median(flux/std)
    print("S/N ratio = ",ston)


    #cont fitting
    intervals = 49
    window = int(speclen/intervals)
    windowwlen = wlen[window]-wlen[0]
    step = 0
    i = 0
    intervalwlen = np.zeros(intervals)
    winpeak = np.zeros(intervals)

    def findmean(array):
        array = np.asarray(array)
        ind = (np.abs(array - np.mean(array))).argmin()
        return ind

    while (step+window) <= speclen:
        windata = flux[step:(step+window)]
        winpeakind = step + findmean(windata)
        winpeak[i] = flux[winpeakind]
        intervalwlen[i] = wlen[winpeakind]
        if i>0 and (intervalwlen[i]-intervalwlen[i-1]) < windowwlen:
            if winpeak[i] <= winpeak[i-1]:
                winpeak[i] = intervalwlen[i] = 0
            else:
                winpeak[i-1] = intervalwlen[i-1] = 0
        step = step + window
        i = i + 1

    winpeak = winpeak[winpeak != 0]
    intervalwlen = intervalwlen[intervalwlen != 0]



    intpol = interpolate.interp1d(intervalwlen, winpeak, kind=1)
    xnew = np.linspace(intervalwlen[0],intervalwlen[-1], num=speclen, endpoint=True)
    ynew = intpol(xnew)

    norm = flux-ynew
    wlim = speclen
    plt.plot(wlen[0:wlim],norm[0:wlim],label=specnames[specind])
    plt.xlabel('wavelength (Angstroms)')
    plt.ylabel('Flux')

plt.legend()


#plotting last processed spectra showing continuum removal
plt.figure()
plt.plot(wlen[0:wlim],flux[0:wlim])
plt.plot(intervalwlen[0:wlim],winpeak[0:wlim],'*')
plt.plot(xnew[0:wlim],ynew[0:wlim],'--')
plt.xlabel('wavelength (Angstroms)')
plt.ylabel('Flux')
plt.show()
