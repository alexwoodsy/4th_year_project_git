import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate
import os

#plt.style.use('mystyle') #path C:\Users\alexw\AppData\Local\Programs\Python\Python37\Lib\site-packages\matplotlib\mpl-data\stylelib

#calculates avgs for inputted array
def findmax(array):
    array = np.asarray(array)
    ind = (np.abs(array - np.max(array))).argmin()
    return ind

def findmed(array):
    array = np.asarray(array)
    ind = (np.abs(array - np.median(array))).argmin()
    return ind

#imports the spectra from the spectra folder
specnames = next(os.walk('Spectra'))[2]
spectot = len(specnames)
#add indexing for spectra in file to allow loop over all

specsample = np.array([0])#indexs of quasars to look at (for later use but added here)

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
    #split the spec in two about lyalpha peak

    forestflux = flux[0:lyalphaind]
    otherflux = flux[lyalphaind:speclen]
    forestwlen = wlen[0:lyalphaind]
    otherwlen = wlen[lyalphaind:speclen]

    forestlen = len(forestflux)
    otherlen = len(otherflux)
    selectlen = np.array([forestlen,otherlen])


    intervalforest, intervalother = 20, 80
    lyalphawidth = 100 # set range around peak for no intervals
    intervals = intervalforest + intervalother

    intervalwlen = np.zeros(intervalforest+intervalother+2)
    winpeak = np.zeros(intervalforest+intervalother+2)

    #loop increments
    windowforest =  int(forestlen/intervalforest)
    windowother =  int(otherlen/intervalother)
    step = windowforest
    i = 0

    while step <= speclen:
        if step <= forestlen:
            window = windowforest
        else:
            window =  windowother

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
        if np.abs(wlen[winpeakind] - wlen[lyalphaind]) > lyalphawidth:
            intervalwlen[i+1] = wlen[winpeakind]
        #elif wlen[winpeakmax] < wlen[lyalphaind] and np.abs(wlen[winpeakind] - wlen[lyalphaind]) < lyalphawidth:
        #    intervalwlen[i+1] = wlen[winpeakind]
        else:
            intervalwlen[i+1] = winpeak[i+1] = 0

        step = step + window
        i = i + 1

#remove zero values made by filtering procedure
    winpeak = winpeak[winpeak != 0]
    intervalwlen = intervalwlen[intervalwlen != 0]

#pad interval with start/end value to allign correctly
    #winpeakmed = step + findmed(windata)
    startind = findmax(flux[0:windowforest])
    winpeak[0],winpeak[-1] = flux[startind],flux[winpeakind]
    intervalwlen[0],intervalwlen[-1] = wlen[0],wlen[-1]

    intpol = interpolate.interp1d(intervalwlen, winpeak, kind=1)
    contfit = intpol(wlen)
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
