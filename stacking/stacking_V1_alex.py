#stacking v1
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

specsample = np.array([0,1000])

gcredhsift = 2
gclyalpha = 1215.67*(1+gcredhsift)
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
