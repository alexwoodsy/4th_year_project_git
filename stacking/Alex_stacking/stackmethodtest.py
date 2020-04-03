import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate
import os

#test data
x1 = np.arange(1,7,0.02)
y1 = np.sin(x1)+1
plt.plot(x1,y1,'.')
x2 = np.arange(1,9,0.05)
y2 = np.sin(x2)
plt.plot(x2,y2,'.')

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
