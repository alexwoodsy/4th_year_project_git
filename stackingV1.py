import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate
import os
import math

plt.style.use('mystyle')

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
# print(specnames)
spectot = len(specnames)
#add indexing for spectra in file to allow loop over all

positions = fits.getdata('PositionsTable.fits',ext=1)
poslen = len(positions)
clusterredshift = np.zeros(poslen)
plate = np.zeros(poslen)
mjd = np.zeros(poslen).astype(str)
fiberid = np.zeros(poslen)
clusternames = np.zeros(poslen).astype(str)
seperation = np.zeros(poslen)

for j in range(0,poslen):
    clusterredshift[j] = positions[j][109]
    clusternames[j] = str(positions[j][105])
    seperation[j] = positions[j][114]
    plate[j] = positions[j][4]
    mjd[j] = str(positions[j][5])
    fiberid[j] = positions[j][6]

# mjdstr = " ".join(str(e) for e in mjd)
# print(clusternames)
mjd = np.array2string(mjd)
print(len(mjd))
# print(type(mjd))

specsample = np.array([0,2,3,4,5]) #indexs of quasars to look at (for later use but added here)

normspecnew = []

for specind in specsample:
    specdirectory = 'Spectra/'+specnames[specind]
    # print(specdirectory)
    data = fits.getdata(specdirectory,ext=1)#import fits image
    speclen = len(data)
    flux = np.zeros(speclen)
    wlen = np.zeros(speclen)

    for i in range(0,speclen):
     flux[i] = data[i][0]
     wlen[i] = 10**(data[i][1])

#meta data extraction to get z:
    fitdata = fits.getdata(specdirectory,ext=2)#import fits image
    # print(fitdata)
    metasize = len(fitdata[0])
    #print(metasize)
    if metasize == 126:
         redshift = fitdata[0][63]
         mjdspec = fitdata[0][56]
    else:
         redshift = fitdata[0][38]
         mjdspec = fitdata[0][29]

    mjdspec = str(mjdspec)
    index = mjd.find(mjdspec)

    clusterz = clusterredshift[index]
    name = clusternames[index]
    sep = seperation[index]
    # print(name)

    lyalphacalc = 1215.67*(1+redshift) #calc lya using redshift
    # lyalphacalc = np.max(flux)
    # print('lyalpha of quasar = ',lyalphacalc)
    lyalphaind = (np.abs(wlen - lyalphacalc)).argmin() #finds index of nearest point in data

    lyalphacluster = 1215.67*(1+clusterz)
    shiftedwlen = wlen - lyalphacluster

#fitting:
    #split the spec in two about lyalpha peak
    forestflux = flux[0:lyalphaind]
    otherflux = flux[lyalphaind:speclen]
    forestwlen = shiftedwlen[0:lyalphaind]
    otherwlen = shiftedwlen[lyalphaind:speclen]
    selecflux = np.array([forestflux,otherflux])
    selecwlen = np.array([forestwlen,otherwlen])

    lyalphawidth = 200 # set range around peak for no intervals

    forestlen = len(forestflux)
    otherlen = len(otherflux)
    selectlen = np.array([forestlen,otherlen])


    intervalforest, intervalother = 40, 12
    intervals = intervalforest + intervalother

    intervalwlen = np.zeros(intervalforest+intervalother+2)
    winpeak = np.zeros(intervalforest+intervalother+2)

    #loop increments
    window =  int(speclen/intervals)
    step = window
    i = 0

    while step <= speclen:
        if step <= forestlen:
            window =  int(forestlen/intervalforest)
        else:
            window =  int(otherlen/intervalother)

        windata = flux[step:(step+window)]
        winpeakmed = step + findmed(windata)
        winpeakmax = step + findmax(windata)
        if shiftedwlen[winpeakmax] < shiftedwlen[lyalphaind]:
            winpeakind = winpeakmax
        elif shiftedwlen[winpeakmax] > shiftedwlen[lyalphaind] and shiftedwlen[winpeakmed] < shiftedwlen[lyalphaind]:
            winpeakind = winpeakmax
        else:
            winpeakind = winpeakmed

        winpeak[i+1] = flux[winpeakind]

#stops slection of interval near lyalpha peak
        if np.abs(shiftedwlen[winpeakind] - shiftedwlen[lyalphaind]) > lyalphawidth:
            intervalwlen[i+1] = shiftedwlen[winpeakind]
        else:
            intervalwlen[i+1] = winpeak[i+1] = 0

        step = step + window
        i = i + 1

#remove zero values made by filtering procedure
    winpeak = winpeak[winpeak != 0]
    intervalwlen = intervalwlen[intervalwlen != 0]

#pad interval with start/end value to allign correctly
    winpeak[0],winpeak[-1] = flux[0],flux[-1]
    intervalwlen[0],intervalwlen[-1] = shiftedwlen[0],shiftedwlen[-1]

    intpol = interpolate.interp1d(intervalwlen, winpeak, kind=1)
    contfit = intpol(shiftedwlen)
    normspec = flux-contfit

    # print(len(normspec))
    normspecnew.append(normspec)
    # normspecnew = normspecnew + normspec
#plotting:
    wlim = speclen
    # plt.plot(shiftedwlen[0:wlim],flux[0:wlim],label=specnames[specind])
    # plt.plot(intervalwlen[0:wlim],winpeak[0:wlim],'*',label='intervals')
    # plt.plot(shiftedwlen[0:wlim],contfit[0:wlim],'--',label='interpolation')
    # plt.plot(shiftedwlen[lyalphaind],flux[lyalphaind],'.',label='lyalpha')
    # plt.xlabel('wavelength (Angstroms)')
    # plt.ylabel('Flux')
    # plt.legend()

    # plt.figure()
    plt.plot(wlen[0:wlim],normspec[0:wlim],label=specnames[specind])
    plt.xlabel('wavelength (Angstroms)')
    plt.ylabel('Flux')
    # plt.legend()

# stacked = np.sum(normspecnew)
# x = zip(normspecnew)
# print(x)
# stacked = sum(x)

# plt.figure()
# plt.plot(shiftedwlen[0:wlim],normspecnew[0:wlim],label=specnames[specind])
# plt.xlabel('wavelength (Angstroms)')
# plt.ylabel('Flux')
# plt.legend()


plt.show()
