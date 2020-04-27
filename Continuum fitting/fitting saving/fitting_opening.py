
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
import os

#plt.style.use('mystyle')

specnames = next(os.walk('Fitted Spectra/'))[2]
spectot = len(specnames)

for spec in specnames:    
    specdirectory = 'Fitted Spectra/' + spec
    data = fits.getdata(specdirectory, ext=1)
    metadata = fits.getdata(specdirectory, ext=1)
    wlen = data.field(0)
    flux = data.field(1)
    medcontfit = data.field(2)
    maxcontfit = data.field(3)


    maxnormspec = flux/maxcontfit

    fig1, ax = plt.subplots(2,1,num='plot')

    ax[0].plot(wlen, flux, label='flux')
    ax[0].plot(wlen, medcontfit, label='med')
    ax[0].plot(wlen, maxcontfit, label='max')
    ax[0].set_xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
    ax[0].set_ylabel(r'$<F>$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')

    ax[1].plot(wlen, maxnormspec, label='normmax')
    ax[1].set_xlabel(r'$\delta$v ($kms^{-1}$)')
    ax[1].set_ylabel(r'MEDIAN $F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')
    ax[1]
    plt.show()











#
