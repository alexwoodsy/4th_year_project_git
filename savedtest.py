import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
import os
#import fitting
import sys
#sys.path.append('C:/Users/jason/GIT/4th_year_project_git/Continuum Fitting')
sys.path.append('C:/Users/alexw/Documents/GitHub/4th_year_project_git/Continuum fitting')


#imports the spectra from the spectra folder
specnames = next(os.walk('Fitted Spectra'))[2]
spectot = len(specnames)

for spec in specnames:
    specdirectory = 'Fitted Spectra/'+spec

    data = fits.getdata(specdirectory,ext=0)
    speclen = len(data)
    fitdata = fits.getdata(specdirectory,ext=1)#read in fits meta data
    metasize = len(fitdata[0])
    #predefine data variables
    flux = np.zeros(speclen)
    wlen = np.zeros(speclen)


    #extract metadata
    for i in range(0,speclen):
        flux[i] = data[i][1]
        wlen[i] = (data[i][0])



    #plotting:
    wlim = speclen
    plt.figure('fitting')
    plt.plot(wlen[0:wlim],flux[0:wlim],label='spec')
    plt.xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
    plt.ylabel(r'$F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')
    plt.show()
