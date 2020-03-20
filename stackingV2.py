#stacking v1
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
import os
#import fitting
import sys
sys.path.append('C:/Users/jason/GIT/4th_year_project_git/Continuum Fitting')
import fittingmethods as fitmeth

#plt.style.use('mystyle')

#imports the spectra from the spectra folder
specnames = next(os.walk('Spectra'))[2]
spectot = len(specnames)

specsample = np.array([0,2,4])

positions = fits.getdata('PositionsTable.fits',ext=1)
poslen = len(positions)
# print(poslen)
# print(positions)
clusterredshift = np.zeros(poslen)
plate = np.zeros(poslen).astype(str)
mjd = np.zeros(poslen).astype(str)
fiberid = np.zeros(poslen).astype(str)
clusternames = np.zeros(poslen).astype(str)
seperation = np.zeros(poslen)

for j in range(0,poslen):
    clusterredshift[j] = positions[j][109]
    clusternames[j] = str(positions[j][105])
    seperation[j] = positions[j][114]
    plate[j] = positions[j][4]
    mjd[j] = str(positions[j][5])
    fiberid[j] = str(positions[j][6])

# mjdstr = " ".join(str(e) for e in mjd)
# print(len(clusternames))
# plate = (plate)
# print(len(plate))
# mjd = np.array2string(mjd)
# print(len(mjd))
# fiberid = np.array2string(fiberid)

normspeckstack = np.zeros(100000)
wlenmin = []
wlenmax = []

for specind in specsample:
    specdirectory = 'Spectra/'+specnames[specind]

    stondata = fits.getdata(specdirectory,ext=1)
    medflux = np.median(stondata[0])
    medivar = np.median(stondata[2])
    std = (np.asarray(medivar))**-0.5
    ston = np.asarray(medflux)/std

    if ston > 0.1:

        fitdata = fits.getdata(specdirectory,ext=2)#import fits image
        # print(fitdata)
        metasize = len(fitdata[0])
        #print(metasize)
        if metasize == 126:
             redshift = fitdata[0][63]
             platespec = fitdata[0][54]
             mjdspec = fitdata[0][56]
             fiberidspec = fitdata[0][57]
        else:
             redshift = fitdata[0][38]
             # platespec = fitdata[0][]
             mjdspec = fitdata[0][29]
             # fiberidspec = fitdata[0][]


        # i for i, x in enumerate(plate) if x==platespec

        # gcredshift = clusterredshift[index2]
        gcredshift = 0
        print(gcredshift)
        gclyalpha = 1215.67*(1+gcredshift)

        wlen, normspec, lyalpha = fitmeth.contfitv6(specind)
        # rfshift = gclyalpha - lyalpha
        # wlenshift = wlen + rfshift
        # wlenshift = wlen - gclyalpha
        wlenshift = wlen/(1+gcredshift)
        wlenmin.append(np.min(wlenshift))
        wlenmax.append(np.max(wlenshift))

        wlenhighres = np.linspace(np.max(wlenmin), np.min(wlenmax), 100000)
        wlenintpol = interpolate.interp1d(wlenshift, normspec, 'linear')
        normspechighres = wlenintpol(wlenhighres)
        normspeckstack = normspeckstack + normspechighres

        # plt.plot(wlenshift,normspec)
        # plt.plot(wlenhighres,normspechighres)

        # print(wlenhighres)

#downsample stacked specind
dsrange = np.linspace(normspeckstack[0], normspeckstack[-1],5000)
dsnormspewcstack = signal.resample(normspeckstack, 5000)
# dsrange = np.linspace(normspeckstack[0], normspeckstack[-1],5000)

#plt.figure()
#plt.plot(wlenhighres, normspeckstack,'.')

plt.figure()
plt.plot(dsrange, dsnormspewcstack)
plt.show()
