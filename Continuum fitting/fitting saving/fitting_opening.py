
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
import os

#plt.style.use('mystyle')

specnames = next(os.walk('Fitted Spectra/'))[2]
spectot = len(specnames)

spec = specnames[0]
specdirectory = 'Fitted Spectra/' + spec
data = fits.getdata(specdirectory, ext=1)
metadata = fits.getdata(specdirectory, ext=1)
wlen = data.field(0)
flux = data.field(1)
medcontfit = data.field(2)
maxcontfit = data.field(3)


maxnormspec = flux/maxcontfit


plt.figure()
plt.plot(wlen, flux, label='flux')
plt.plot(wlen, medcontfit, label='med')
plt.plot(wlen, maxcontfit, label='max')
plt.plot(wlen, maxnormspec, label='normmax')
plt.show()











#
