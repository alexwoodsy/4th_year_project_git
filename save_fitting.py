import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
import os
#import fitting
import sys
sys.path.append('C:/Users/jason/GIT/4th_year_project_git/Continuum Fitting')
import fitting_jason as fitmeth

specnames = next(os.walk('Spectra'))[2]
# zlim = gcredshift+0.15
stonlim = 1

for spec in specnames:
    specdirectory = 'Spectra/'+spec
    spec = [spec]
    wlen, normspec, wlenlineind, redshift,stackstatus = fitmeth.contfitv7(spec, stonlim, showplot = False, showerror = False)

    # out = open("fitted",'a')
    # print(specdirectory, wlen, normspec, wlenlineind, redshift,file=out)
    # out.close()
    # print(spec)
    # fits.writeto(spec, wlen, header=wlen.header_model, clobber=False)
    # fits.writeto(spec, normspec, header=normspec.header_model, clobber=False)
    # fits.writeto(spec, wlenlineind, header=wlenlineind.header_model, clobber=False)
    # fits.writeto(spec, redshift, header=redshift.header_model, clobber=False)

    c1 = fits.Column(name='Wavelength', array=wlen, format='K')
    c2 = fits.Column(name='Flux', array=normspec, format='K')
    c3 = fits.Column(name='Wavelength Index', array=wlenlineind, format='K')
    c4 = fits.Column(name='Spectra Redshift', array=redshift, format='K')
    t = fits.BinTableHDU.from_columns([c1, c2, c3])
    t.writeto('prefitted.fits')
