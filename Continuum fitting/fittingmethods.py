#all methods collated to be called as continuum fittind function
#without plotting and contiuum reading in

import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate
import os

#imports the spectra from the spectra folder
specnames = next(os.walk('Spectra'))[2]
spectot = len(specnames)



####------V6 multi-interval+model method------#####

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


def contfitv4(specind):
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


    intervalforest, intervalother = 20, 70
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

    return wlen, normspec

def contfitv5(specind):
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
    #print('lyalpha = ',lyalphacalc)
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


    intervalforest = 20
    # set range around peak for no intervals
    lyalphalimind = 25
    forestlen = forestlen-lyalphalimind


    intervalwlen = np.zeros(intervalforest+1)
    winpeak = np.zeros(intervalforest+1)

    #loop increments
    window =  int(forestlen/intervalforest)
    step = window
    i = 0

    while step <= forestlen:
        windata = flux[step:(step+window)]
        winpeakmax = step + findmax(windata)
        winpeakind = winpeakmax
        winpeak[i+1] = flux[winpeakind]
#stops slection of interval near lyalpha peak
        intervalwlen[i+1] = wlen[winpeakind]
        step = step + window
        i = i + 1

#remove zero values made by filtering procedure
    winpeak = winpeak[winpeak != 0]
    intervalwlen = intervalwlen[intervalwlen != 0]

#pad interval with start/end value to allign correctly
    #winpeakmed = step + findmed(windata)
    startind = findmax(flux[0:window])
    winpeak[0] = flux[startind]
    intervalwlen[0] = wlen[0]

    intpol = interpolate.interp1d(intervalwlen, winpeak, kind=1)
    interpolfit = intpol(forestwlen[0:forestlen])
    contfit = model
    contfit[0:forestlen] = interpolfit


    normspec = flux-contfit
    return wlen, normspec


def contfitv6(specind):
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
    #print('lyalpha = ',lyalphacalc)
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
    startind = findmax(flux[0:window])
    winpeak[0] = flux[startind]
    intervalwlen[0] = wlen[0]

    intpol = interpolate.interp1d(intervalwlen, winpeak, kind=1)
    interpolfit = intpol(forestwlen)
    contfit = model
    contfit[0:forestlen] = interpolfit

    normspec = flux-contfit
    return wlen, normspec
