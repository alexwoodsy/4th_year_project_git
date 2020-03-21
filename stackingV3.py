#stacking v2
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
import os
#import fitting
import sys
sys.path.append('C:/Users/jason/GIT/4th_year_project_git/Continuum Fitting')

# sys.path.append('C:/Users/alexw/Documents/GitHub/4th_year_project_git/Continuum fitting')
#path for other pc sys.path.append('C:/Users/alexw/OneDrive/Documents/University work/4th year work/Main project/4th_year_project_git/Continuum fitting')
import fitting_v7 as fitmeth

# plt.style.use('mystyle')

#imports the spectra from the spectra folder
specnames = next(os.walk('Spectra'))[2]
spectot = len(specnames)

#set up sample to be looped over
#read in data linking cluster to qso's
posdata = fits.getdata('PositionsTable.fits',ext=1)
poslen = len(posdata)
clusterredshift = np.zeros(poslen)
clusternames = np.zeros(poslen).astype(str)
specfilename = np.zeros(poslen).astype(str)

normspeckstack = np.zeros(100000)
wlenmin = []
wlenmax = []

for j in range(0,poslen):
    #get cluster info
    clusterredshift[j] = posdata[j][109]
    clusternames[j] = str(posdata[j][105])
    #get assocaited id to find spec in folder
    plate = str(posdata[j][4]).zfill(4)
    mjd = str(posdata[j][5]).zfill(5)
    fiberid = str(posdata[j][6]).zfill(4)
    specfilename[j] = 'spec-'+plate+'-'+mjd+'-'+fiberid+'.fits'

#read in carla
carladata = fits.getdata('CARLA/CARLA_table_Aug_2013.fits',ext=1)
carlalen = len(carladata)
carlanames = np.zeros(carlalen).astype(str)
for i in range(0,carlalen):
    carlanames[i] = str(carladata[i][0])

#select carla agn that were selected in filtering (in pos table) and associated spec index in table
match =  []

for i in carlanames:
    c = 0
    for k in clusternames:
        if i == k:
            c = c + 1
    if c > 1:
        match.append(i)

#get spec names for a carla target - do stacking

for carlaselect in match[0:1]:
    specmatch = []
    for i in range(0,poslen):
        if clusternames[i] == carlaselect and specnames[i][-4:] == 'fits':
            specmatch.append(specnames[i])

#stack here - specmatch = qso to stack
    wlen, normspec, lyalpha = fitmeth.contfitv7(specmatch[0:8], showplot = True)

    clusterredshift = 0
    wlenshift = wlen/(1+clusterredshift)
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
plt.xlabel('wavelength (Angstroms)')
plt.ylabel('Flux')
plt.title('Stacked?')
plt.show()

#
