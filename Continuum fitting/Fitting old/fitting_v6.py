#top x% method for fitting continuum in the forest
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate
import os

#plt.style.use('mystyle') #path C:\Users\alexw\AppData\Local\Programs\Python\Python37\Lib\site-packages\matplotlib\mpl-data\stylelib

#calculates points slected in interval
def findmax(array):
    array = np.asarray(array)
    ind = (np.abs(array - np.max(array))).argmin()
    return ind

def findmed(array):
    array = np.asarray(array)
    ind = (np.abs(array - np.median(array))).argmin()
    return ind

def findpct(array,pct):
    sortedarray = np.argsort(array)
    selectedvals = sortedarray[-pct:]
    return selectedvals


# imports the spectra from the spectra folder
specnames = next(os.walk('Spectra'))[2]
spectot = len(specnames)
# add indexing for spectra in file to allow loop over all

specsample = np.array([1000])#indexs of quasars to look at (for later use but added here)

for specind in specsample:
    specdirectory = 'Spectra/'+specnames[specind]
    #print(specdirectory)
    data = fits.getdata(specdirectory,ext=1)#import fits image
    speclen = len(data)
    flux = np.zeros(speclen)
    wlen = np.zeros(speclen)
    model = np.zeros(speclen)

    for i in range(0,speclen):
     flux[i] = data[i][0]
     wlen[i] = 10**(data[i][1])
     model[i] = data[i][7]


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
    #split the spec in two about lyalpha peak

    forestflux = flux[0:lyalphaind]
    otherflux = flux[lyalphaind:speclen]
    forestwlen = wlen[0:lyalphaind]
    otherwlen = wlen[lyalphaind:speclen]

    forestlen = len(forestflux)
    otherlen = len(otherflux)
    selectlen = np.array([forestlen,otherlen])

    intervalwlen = np.array([0])
    winpeak = np.array([0])
    pct = 15
    #loop increments
    winnum = 10
    window = int(forestlen/winnum)
    step = 0

    while step <= forestlen+winnum:
        windata = flux[step:(step+window)]
        winpeakind = step + findpct(windata,pct)
        winpeak = np.append(winpeak,flux[winpeakind])
        intervalwlen = np.append(intervalwlen,wlen[winpeakind])
        step = step + window


#pad interval with start/end value to allign correctly
    #winpeakmed = step + findmed(windata)
    startind = findmed(flux[0:window])
    winpeak[0] = flux[startind]
    intervalwlen[0] = wlen[0]


    intpol = interpolate.interp1d(intervalwlen, winpeak, kind=1)
    interpolfit = intpol(forestwlen)
    contfit = model
    contfit[0:forestlen] = interpolfit


    normspec = flux-contfit


#plotting:
    wlim = speclen
    plt.figure(1)
    #plt.title('continuum fitting')
    plt.plot(wlen[0:wlim],flux[0:wlim],label=specnames[specind][:20])
    plt.plot(intervalwlen[0:wlim],winpeak[0:wlim],'*',label='intervals')
    plt.plot(wlen[0:wlim],contfit[0:wlim],'--',label='interpolation')
    plt.plot(wlen[lyalphaind],flux[lyalphaind],'.',label=r'Ly$\alpha$')
    plt.xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
    plt.ylabel(r'$F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')
    if specind == specsample[0]:
        plt.legend()

    plt.figure(2)
    #plt.title('continuum removed')
    plt.plot(wlen[0:wlim],normspec[0:wlim],label=specnames[specind][:20])
    plt.xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
    plt.ylabel(r'$F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')
    plt.legend()

plt.show()
