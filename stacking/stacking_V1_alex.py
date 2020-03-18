#stacking v1
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate
import os
#import fitting
import sys
sys.path.append('C:/Users/alexw/Documents/GitHub/4th_year_project_git/Continuum fitting')
import fittingmethods as fitmeth

#plt.style.use('mystyle') #path C:\Users\alexw\AppData\Local\Programs\Python\Python37\Lib\site-packages\matplotlib\mpl-data\stylelib

#imports the spectra from the spectra folder
specnames = next(os.walk('Spectra'))[2]
spectot = len(specnames)

specsample = np.array([0,2,1000])


#test data
x1 = np.arange(1,7,0.02)
y1 = np.sin(x1)
plt.plot(x1,y1)
x2 = np.arange(1,9,0.05)
y2 = np.sin(x2)
plt.plot(x2,y2)

shift = 0 - np.max(y1)
xshift = x1 + shift
xhighres = np.linspace(np.min(xshift), np.max(xshift), 1000)
xintpol = interpolate.interp1d(xshift, y1, 'linear')
y1highres = xintpol(xhighres)

shift = 0 - np.max(y2)
xshift = x2 + shift
xhighres = np.linspace(np.min(xshift), np.max(xshift), 1000)
xintpol = interpolate.interp1d(xshift, y2, 'linear')
y2highres = xintpol(xhighres)

yout = y1highres-y2highres
plt.plot(xhighres,yout)
plt.show()




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
plt.figure()
plt.plot(wlenhighres, normspeckstack)
plt.show()
