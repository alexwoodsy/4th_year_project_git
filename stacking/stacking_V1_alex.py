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
i = 0
for specind in specsample:
    wlen, normspec, lyalpha = fitmeth.contfitv6(specind)
    i = i + 1
    print(i)



    #wlim = len(normspec)
    #plt.plot(wlen[0:wlim],normspec[0:wlim],label=specnames[specind][:20])
    #plt.xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
    #plt.ylabel(r'$F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')
    #plt.legend()
#plt.show()
