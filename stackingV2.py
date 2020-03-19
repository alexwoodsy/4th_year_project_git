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

specsample = np.array([0,2,1000])

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

normspeckstack = np.zeros(100000)

for specind in specsample:
    specdirectory = 'Spectra/'+specnames[specind]
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
    gcredshift = clusterredshift[index]
    gclyalpha = 1215.67*(1+gcredshift)

    wlen, normspec, lyalpha = fitmeth.contfitv6(specind)
    rfshift = gclyalpha - lyalpha
    wlenshift = wlen + rfshift
    # wlenshift = wlen - gclyalpha
    wlenhighres = np.linspace(np.min(wlenshift), np.max(wlenshift), 100000)
    wlenintpol = interpolate.interp1d(wlenshift, normspec, 'linear')
    normspechighres = wlenintpol(wlenhighres)
    normspeckstack = normspeckstack + normspechighres


#downsample stacked specind
dsrange = np.linspace(normspeckstack[0], normspeckstack[-1],5000)
dsnormspewcstack = signal.resample(normspeckstack, 5000)
dsrange = np.linspace(normspeckstack[0], normspeckstack[-1],5000)

#plt.figure()
#plt.plot(wlenhighres, normspeckstack,'.')

plt.figure()
plt.plot(dsrange, dsnormspewcstack)
plt.show()
