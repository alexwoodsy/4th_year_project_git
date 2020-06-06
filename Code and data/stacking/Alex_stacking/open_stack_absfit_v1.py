#stacking v2
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
from scipy.optimize import curve_fit as cf
import os, time
#import fitting
import sys
sys.path.append('C:/Users/alexw/Documents/GitHub/4th_year_project_git/Continuum fitting')
#path for other pc sys.path.append('C:/Users/alexw/OneDrive/Documents/University work/4th year work/Main project/4th_year_project_git/Continuum fitting')
import fitting_v9 as fitmeth
#plt.style.use('mystyle') #path C:\Users\alexw\AppData\Local\Programs\Python\Python37\Lib\site-packages\matplotlib\mpl-data\stylelib

#open stacking data
stackpath = 'stacking/figures/Stacking data/'
run_name = 'rad_bins_(3000_to_3500).fits'

hdul = fits.open(stackpath+run_name)  # open a FITS file
header = hdul[1].header
data = hdul[1].data
hdul.close()
colnumber = header['TFIELDS']
rownumber = header['NAXIS2']

wlen = np.zeros(rownumber)
meanflux = np.zeros(rownumber)
medflux = np.zeros(rownumber)
wlen[0:] = data.field(0)
meanflux[0:] = data.field(1)
medflux[0:] = data.field(2)

###absorbtion line fitting####
#fitting needs initial data so extract data about abs line
def guassian(x, amp, mean, std):
    return amp*np.exp(-((x-mean)**2)/(2*(std)**2))

def findval(array,val):
    array = np.asarray(array)
    ind = np.abs(array - val).argmin()
    return ind

fig0, ax = plt.subplots(2,1,num=run_name+'Absorption line fitting')

abslineind = findval(wlen, 1215.67)
datarange = np.arange(abslineind - 12, abslineind + 15)
plotrange = np.arange(abslineind - 50, abslineind + 50)
absflux = meanflux[datarange]
abswlen = wlen[datarange]

popt, pcov = cf(guassian, abswlen, absflux, bounds =([-np.inf,1215.67-0.5,-np.inf],[np.inf,1215.67+0.5,np.inf]))
meanamp, meanmean, meanstd = popt
ax[0].plot(wlen[plotrange], meanflux[plotrange])
ax[0].plot(abswlen, guassian(abswlen, *popt), 'r-',label='fitting parmaters: amp=%5.3f, mean=%5.3f, std=%5.3f' % tuple(popt))
ax[0].set_xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
ax[0].set_ylabel(r'$<F>$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')
ax[0].legend()

absflux = medflux[datarange]
abswlen = wlen[datarange]
popt, pcov = cf(guassian, abswlen, absflux, bounds =([-np.inf,1215.67-0.5,-np.inf],[np.inf,1215.67+0.5,np.inf]))
medamp,medmean,medstd = popt
ax[1].plot(wlen[plotrange], medflux[plotrange])
ax[1].plot(abswlen, guassian(abswlen, *popt), 'r-',label='fitting parmaters: amp=%5.3f, mean=%5.3f, std=%5.3f' % tuple(popt))
ax[1].set_xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
ax[1].set_ylabel(r'MEDIAN $F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')
ax[1].legend()

#plot relative vel graph and find delv of line from rest gc ly alpha
fig1, ax = plt.subplots(2,1,num=run_name+'velocity Absorption line plot')
c = 299792.458
lam = wlen
lam_em = 1215.67
vrel = c*((lam  - lam_em)/lam_em)

abslineind = findval(vrel, 0)
plotrange = np.arange(abslineind - 50, abslineind + 50)

ax[0].plot(vrel, meanflux)
ax[0].set_xlabel(r'$\delta$v ($kms^{-1})$')
ax[0].set_ylabel(r'$<F>$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')

ax[1].plot(vrel, medflux)
ax[1].set_xlabel(r'$\delta$v ($kms^{-1}$)')
ax[1].set_ylabel(r'MEDIAN $F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')


fig3, ax = plt.subplots(2,1,num=run_name+' stacked clusters')
ax[0].plot(wlen[:-300], meanflux[:-300])
ax[0].plot(np.array([(1215.67),(1215.67)]),np.array([np.min(meanflux),np.max(meanflux)]),'--')
ax[0].text(1215.67, np.min(meanflux), r' Ly$\alpha$ absorbtion for multi-carla stack')
ax[0].set_xlabel(r' $\lambda$ ($\mathrm{\AA}$)')
ax[0].set_ylabel(r'$<F>$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')

ax[1].plot(wlen[:-300], medflux[:-300])
ax[1].plot(np.array([(1215.67),(1215.67)]),np.array([np.min(medflux),np.max(medflux)]),'--')
ax[1].text(1215.67, np.min(medflux), r' Ly$\alpha$ absorbtion for multi-carla stack')
ax[1].set_xlabel(r' $\lambda$ ($\mathrm{\AA}$)')
ax[1].set_ylabel(r'MEDIAN $F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')



plt.show()
