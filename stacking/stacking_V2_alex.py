#stacking v2
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
import os
#import fitting
import sys
sys.path.append('C:/Users/alexw/Documents/GitHub/4th_year_project_git/Continuum fitting')
#path for other pc sys.path.append('C:/Users/alexw/OneDrive/Documents/University work/4th year work/Main project/4th_year_project_git/Continuum fitting')
import fittingmethods as fitmeth

#plt.style.use('mystyle') #path C:\Users\alexw\AppData\Local\Programs\Python\Python37\Lib\site-packages\matplotlib\mpl-data\stylelib

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
for k in range(0,carlalen):
    c = 0
    for namecheck in clusternames:
        if namecheck == carlanames[k]:
            c = c + 1
    if c > 0:
        match.append(namecheck)

#get spec names for a carla target
carlatarget = match[0]
matchspecname =  []
for i in range(0,poslen):
    if clusternames[i] == carlatarget:
        matchspecname.append(specfilename[i])
        gcredhsift = clusterredshift[i]




gcredshift = 2
gclyalpha = 1215.67*(1+gcredhsift)
specsample = np.array(])


#loops overspectra, upsamling - alligning stacking and then downsampling
normspeckstack = np.zeros(100000)
for specind in specsample:
    wlen, normspec, lyalpha = fitmeth.contfitv6(specind)
    rfshift = gclyalpha - lyalpha
    wlenshift = wlen + rfshift
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
