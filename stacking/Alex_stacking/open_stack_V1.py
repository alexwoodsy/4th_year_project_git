#stacking v2
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
import os, time
#import fitting
import sys
sys.path.append('C:/Users/alexw/Documents/GitHub/4th_year_project_git/Continuum fitting')
#path for other pc sys.path.append('C:/Users/alexw/OneDrive/Documents/University work/4th year work/Main project/4th_year_project_git/Continuum fitting')
import fitting_v9 as fitmeth
plt.style.use('mystyle') #path C:\Users\alexw\AppData\Local\Programs\Python\Python37\Lib\site-packages\matplotlib\mpl-data\stylelib

#open stacking data
stackpath = 'stacking/figures/Stacking data/'
run_name = 'testsave.fits'

hdul = fits.open(stackpath+run_name)  # open a FITS file
header = hdul[1].header
data = hdul[1].data
hdul.close()
colnumber = header['TFIELDS']
rownumber = header['NAXIS2']

wlen = np.zeros([rownumber,colnumber])
flux = np.zeros([rownumber,colnumber])
gcnames = []
for c in range(0,colnumber,2):
    gcname = header['TTYPE'+str(c+1)]
    gcnames.append(gcname[:-5])
    wlen[0:, c] = data.field(c)
    flux[0:, c+1] = data.field(c+1)
    plt.plot(wlen[0:, c],flux[0:, c+1],label=gcname)



plt.show()




#cut down zero padding
# gcstart = (np.abs(wlenmultistack - gcwlenmin)).argmin()
# gcend = (np.abs(wlenmultistack - gcwlenmax)).argmin()
# wlenmultistack = wlenmultistack[gcstart:gcend]
# multistack = multistack[gcstart:gcend]

#plotting and data manipulation for output

# meanstack = multistack/specstacktot
# medianstack = multistack
# #
#
# #plot line in unstacked graph
# plt.figure('multi-carla stack + uncombined')
# plt.plot(np.array([(1215.67),(1215.67)]),np.array([150,-150]),'--')
# plt.text(1215.67, 150, r' Ly$\alpha$ absorbtion for '+ str(carlatot) +' carla with '+ str(specstacktot) +' spectra')
#
# plt.figure('stacked clusters')
# plt.subplot(2,1,1)
# plt.plot(wlenmultistack[:-300], meanstack[:-300])
# plt.plot(np.array([(1215.67),(1215.67)]),np.array([np.min(meanstack),np.max(meanstack)]),'--')
# plt.text(1215.67, np.min(meanstack), r' Ly$\alpha$ absorbtion for multi-carla stack with '+ str(specstacktot) +' spectra')
# plt.xlabel(r' $\lambda$ ($\mathrm{\AA}$)')
# plt.ylabel(r'$<F>$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')
#
#
# plt.subplot(2,1,2)
# plt.plot(wlenmultistack[:-300], medianstack[:-300])
# plt.plot(np.array([(1215.67),(1215.67)]),np.array([np.min(meanstack),np.max(meanstack)]),'--')
# plt.text(1215.67, np.min(meanstack), r' Ly$\alpha$ absorbtion for multi-carla stack with '+ str(specstacktot) +' spectra')
# plt.xlabel(r' $\lambda$ ($\mathrm{\AA}$)')
# plt.ylabel(r'$<F>$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')
#
# plt.show()
