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
    if spec[-4:] == 'fits':
        specdirectory = 'Spectra/' + spec
        spec = [spec]
        wlen, normspec, redshift, stackstatus = fitmeth.contfitv7(spec, stonlim, showplot = False, showerror = False)
        spec = str(spec)


        c1 = fits.Column(name='Wavelength', array=wlen, format='K')
        c2 = fits.Column(name='Flux', array=normspec, format='K')
        c4 = fits.Column(name='Spectra Redshift', array=redshift, format='K')
        t = fits.BinTableHDU.from_columns([c1, c2])
        # different headers
        meta = fits.BinTableHDU.from_columns(c4)
        # forest s/n save
        hdul = fits.HDUList([t,meta])
        #
        hdul.writeto('Fitted Spectra/' + spec[2:22] + '-prefitted.fits',overwrite = True)
        # fits.writeto(hdul,'Fitted Spectra/' + spec[2:22] + '-prefitted.fits',overwrite = True)
