import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
from scipy.optimize import curve_fit as cf
import os, random
#import fitting
import sys

plt.style.use('mystyle') #path C:\Users\alexw\AppData\Local\Programs\Python\Python37\Lib\site-packages\matplotlib\mpl-data\stylelib

def findval(array,val):
    array = np.asarray(array)
    ind = np.abs(array - val).argmin()
    return ind



#gc = J092058.46+444154.0
#'spec-0832-52312-0566-prefitted.fits'
#'spec-0832-52312-0113-prefitted.fits'
#'spec-0833-52314-0194-prefitted.fits'

specstacknames = ['spec-0832-52312-0566-prefitted.fits',
'spec-0832-52312-0113-prefitted.fits','spec-0833-52314-0194-prefitted.fits']

######-----stacking code----#######
fillval = np.nan
runlen = len(specstacknames)
highreslen = 50000
cutinds = []
normspecstore = np.empty([runlen,highreslen])
contspecstore = np.empty([runlen,highreslen])
normspecstore[:] = fillval
contspecstore[:] = fillval

wlenhighres = np.linspace(500, 4500, highreslen)
wlenmin = 10000
wlenmax = 0
stackstatus = []
specnumber = 0



for spec in specstacknames:
    specdata = fits.getdata('Fitted Spectra/'+spec,ext=1)#import fits image
    wlen = specdata.field(0)
    flux = specdata.field(1)
    contfit = specdata.field(6)
    specmetadata = fits.getdata('Fitted Spectra/'+spec,ext=2)#import fits image
    redshift = specmetadata[0][0]
    gcname = str(specmetadata[0][5])
    gcredshift = specmetadata[0][6]
    print(gcredshift)

    cutstart = findval(wlen,(1030*(1+gcredshift)))
    cutend = findval(wlen,(1600*(1+gcredshift)))
    if cutend != 0 and cutstart != 0:
        cut = np.arange(cutstart,cutend)
        contfit = contfit[cut]
        flux = flux[cut]
        wlen = wlen[cut]

    wlenshift = wlen/(1+gcredshift) #+np.random.random(1)-0.5
    normspec = flux/contfit

    # plt.plot(wlenshift,flux)
    # plt.plot(wlenshift,contfit)
    #plt.plot(wlenshift,normspec)
    #plt.show()

    wlenintpol = interpolate.interp1d(wlenshift, normspec, 'linear', bounds_error=False, fill_value=fillval)
    contintpol = interpolate.interp1d(wlenshift, contfit, 'linear', bounds_error=False, fill_value=fillval)
    if wlenshift[0] < wlenmin:
        wlenmin = wlenshift[0]
    if wlenshift[-1] > wlenmax:
        wlenmax = wlenshift[-1]




    normspecstore[specnumber, 0:] = wlenintpol(wlenhighres)
    contspecstore[specnumber, 0:] = contintpol(wlenhighres)

    specnumber = specnumber + 1


#stacking data only if there are spec to stack

#cut extra zeropadding
start = (np.abs(wlenhighres - wlenmin)).argmin()+10 #accounts for inaccuracy in cut
end = (np.abs(wlenhighres - wlenmax)).argmin() -10
wlenhighres = wlenhighres[start:end]
#cut downcols / remove empty rows
#normspecstore = np.delete(normspecstore, cutinds, axis = 0)
normspecstore = normspecstore[0:, start:end]
meanspec = np.nanmean(normspecstore, axis=0)
medspec = np.nanmedian(normspecstore, axis=0)
print(normspecstore)



fig = plt.figure('stacking_example')
# set height ratios for sublots
gs = plt.GridSpec(3, 1, height_ratios=[3, 2, 2])
# the fisrt subplot
ax0 = plt.subplot(gs[0])
#ax0.set_yscale("log")
ax0.plot(wlenhighres,normspecstore[0, 0:],alpha=0.6,color='#ff0000')
ax0.plot(wlenhighres,normspecstore[1, 0:],alpha=0.6,color='#00ff00')
ax0.plot(wlenhighres,normspecstore[2, 0:],alpha=0.6,color='#0000ff')

ax1 = plt.subplot(gs[1], sharex = ax0)
ax1.plot(wlenhighres, meanspec,color='#000000')
#ax1.set_ylabel(r'$\tilde{F}$')

ax2 = plt.subplot(gs[2], sharex = ax0, sharey = ax1)
ax2.plot(wlenhighres, medspec,color='#808080')
ax2.set_xlabel(r'$\lambda$ ($\mathrm{\AA}$)')

# remove vertical gap between subplots
plt.subplots_adjust(hspace=.0)
#labels
# Set common labels
#fig.text(0.5, 0.03, r'$\lambda$ ($\mathrm{\AA}$)', ha='center', va='center')
fig.text(0.05, 0.5, r'$\tilde{F}$', ha='center', va='center', rotation='vertical')
#fig.text(0.05, 0.5, r'$F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$', ha='center', va='center', rotation='vertical')
plt.show()
