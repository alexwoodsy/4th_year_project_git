#metadata extraction of redshift and s to n
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
import os

specnames = next(os.walk('Spectra'))[2] #dir is your directory path as string
spectot = len(specnames)

#add indexing for epctra in file to allow loop over all
number = 10
redshift = np.zeros(number)
snmedian = np.zeros(number)

for i in range(0,number):
    specdirectory = 'Spectra/'+specnames[i]
    print(specdirectory)

    fitdata = fits.getdata(specdirectory,ext=2)#import fits image

    metasize = len(fitdata[0])
    #print(metasize)
    if metasize == 126:
        redshift[i] = fitdata[0][63]
        snmedian[i] = np.max(fitdata[0][84])
        #for some reason they take the value from the 3rd filter
    else:
        redshift[i] = fitdata[0][38]
        snmedian[i] = np.median(fitdata[0][58])

print(snmedian)
