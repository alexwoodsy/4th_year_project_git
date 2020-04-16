import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate
import os

plt.style.use('mystylehistogram')


specdirectory = 'metatable.fits'
fitdata = fits.getdata(specdirectory,ext=1)#read in fits meta data
fitlen = len(fitdata)

specfilename = np.zeros(fitlen).astype(str)
redshift = np.zeros(fitlen)
stonall = np.zeros(fitlen)
stonforest = np.zeros(fitlen)
lyalpha = np.zeros(fitlen)
lyalphaind = np.zeros(fitlen)
gcname = np.zeros(fitlen).astype(str)
gcredshift = np.zeros(fitlen)
gc_qso_sep = np.zeros(fitlen)
gclyalpha = np.zeros(fitlen)
gclyalphaind = np.zeros(fitlen)
stackmsg = np.zeros(fitlen).astype(str)

for i in range(0,fitlen):
    #extract metadata
    specfilename[i] = fitdata[i][0]
    redshift[i]= fitdata[i][1]
    stonall[i] = fitdata[i][2]
    stonforest[i] = fitdata[i][3]
    lyalpha[i] = fitdata[i][4]
    lyalphaind[i] = fitdata[i][5]
    gcname[i] = fitdata[i][6]
    gcredshift[i] = fitdata[i][7]
    gc_qso_sep[i] = fitdata[i][8]
    gclyalpha[i] = fitdata[i][9]
    gclyalphaind[i] = fitdata[i][10]
    stackmsg[i] = fitdata[i][11]

ind = (np.argmax(redshift))

plt.figure('s/n hist')
xvar = [stonall,stonforest]
step = 1
binning = np.arange(int(np.min(xvar)), int(np.max(xvar)+step) , step)
labelnames = [r'Spectrum Range', r'Ly$\alpha$ forest']
n, bins, patches = plt.hist(x = xvar, bins = binning, log = True, alpha=1, rwidth = 1, histtype = 'bar', label = labelnames)

#plt.grid(axis='y', alpha=0.4)
plt.xlabel(r'$(S/N)$')
plt.xticks(binning)
plt.ylabel(r'$N_{QSO}$')
plt.legend()

plt.figure('rad hist')
xvar = gc_qso_sep
step = 50
print(np.min(xvar))
print(np.max(xvar))
binning = np.arange(int(np.min(xvar)), int(np.max(xvar)+step) , step)
labelnames = [r'$R_{\perp}$']
n, bins, patches = plt.hist(x = xvar, bins = binning, log = True, alpha=1, rwidth = 0.9, histtype = 'bar', label = labelnames)

#plt.grid(axis='y', alpha=0.4)
plt.xlabel(r'$R_{\perp}$ (arcseconds)')
tickrange = np.arange(int(np.min(xvar)), int(np.max(xvar)+step) , step*4)
plt.xticks(tickrange)
plt.ylabel(r'$N_{QSO}$')
plt.legend()



plt.show()














#
