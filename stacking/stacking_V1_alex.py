#stacking v1
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate
import os
#import fitting
import sys
sys.path.append('C:/Users/alexw/OneDrive/Documents/University work/4th year work/Main project/4th_year_project_git/Continuum fitting')
import fittingmethods as fitmeth

#plt.style.use('mystyle') #path C:\Users\alexw\AppData\Local\Programs\Python\Python37\Lib\site-packages\matplotlib\mpl-data\stylelib

#imports the spectra from the spectra folder
specnames = next(os.walk('Spectra'))[2]
spectot = len(specnames)
specsample = np.array([0,2,1000])



fitsloc = 'stacking/contfitraw.fits'

col = fits.Column(name=specnames[0], array=np.array([1, 5]), format='K')
tablehdu = fits.BinTableHDU.from_columns([col])

wlen = np.array([1,2,3,4,5,6])
normspec = np.array([1,2,3,4,5,6])

for specind in specsample[1:]:
    data = wlen
    col = fits.Column(name=specnames[specind], array=data, format='K')
    tablehdunew = fits.BinTableHDU.from_columns([col])
    new_columns = tablehdu.columns + tablehdunew.columns
    tablehdu = fits.BinTableHDU.from_columns(new_columns)





tablehdu.writeto(fitsloc, clobber = True)# initially makes file (not needed after 1st run)
rawdata = fits.getdata(fitsloc,ext=1)
print(rawdata.field(0))
print(rawdata.field(1))
print(rawdata.field(2))



    #wlim = len(normspec)
    #plt.plot(wlen[0:wlim],normspec[0:wlim],label=specnames[specind][:20])
    #plt.xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
    #plt.ylabel(r'$F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')
    #plt.legend()
#plt.show()
