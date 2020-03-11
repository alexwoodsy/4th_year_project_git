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
    #model = np.zeros(speclen)
    #ivar = np.zeros(speclen)

    for i in range(0,speclen):
     flux[i] = data[i][0]
     #ivar[i] = data[i][2]
     #model[i] = data[i][7]
     wlen[i] = 10**(data[i][1])

 #meta data extraction to get z:
    fitdata = fits.getdata(specdirectory,ext=2)#import fits image
    metasize = len(fitdata[0])
    #print(metasize)
    if metasize == 126:
         redshift = fitdata[0][63]
    else:
         redshift = fitdata[0][38]


    lyalphacalc = 1026.67*(1+redshift)#calc lya using redshift
    print('lyalpha = ',lyalphacalc)
    lyalphaind = (np.abs(wlen - lyalphacalc)).argmin()#finds index of nearest point in data


#fitting: does fit to max of ly alpha and then median for rest of spec
    intervals = 60

    forestflux = flux[0:lyalphaind]
    otherflux = flux[lyalphaind:speclen]
    forestwlen = wlen[0:lyalphaind]
    otherwlen = wlen[lyalphaind:speclen]

    window = int(speclen/intervals)

    windowwlen = wlen[window]-wlen[0]
    print(windowwlen)
    step = 0
    i = 0
    intervalwlen = np.zeros(intervals)
    winpeak = np.zeros(intervals)

    def findmax(array):
        array = np.asarray(array)
        ind = (np.abs(array - np.max(array))).argmin()
        return ind

    def findmed(array):
        array = np.asarray(array)
        ind = (np.abs(array - np.med(array))).argmin()
        return ind

    while (step+window) <= speclen:
        windata = flux[step:(step+window)]
        winpeakind = step + findmax(windata)
        winpeak[i] = flux[winpeakind]
        intervalwlen[i] = wlen[winpeakind]
        step = step + window
        i = i + 1

    intpol = interpolate.interp1d(intervalwlen, winpeak, kind=1)
    xnew = np.linspace(intervalwlen[0],intervalwlen[-1], num=speclen, endpoint=True)
    ynew = intpol(xnew)

#plotting:
    wlim = speclen
    plt.plot(wlen[0:wlim],flux[0:wlim],label=specnames[specind])
    plt.plot(wlen[0:wlim],ynew[0:wlim],label=specnames[specind])
    #plt.plot(intervalwlen[0:wlim],winpeak[0:wlim],'*',label='intervals')
    #plt.plot(wlen[0:wlim],ynew[0:wlim],'--',label='interpolation')
    #plt.plot(wlen[lyalphaind],flux[lyalphaind],'.',label='lyalpha')
    plt.xlabel('wavelength (Angstroms)')
    plt.ylabel('Flux')


plt.legend()
plt.show()
