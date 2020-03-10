#trgt_smple_slct
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os

specnames = next(os.walk('Spectra'))[2] #dir is your directory path as string
spectot = len(specnames)

#sample selector
specsample = np.array([0,1,2])
sslen = len(specsample)



#add indexing for epctra in file to allow loop over all

redshift = np.zeros(sslen)
snmedian = np.zeros(sslen)

for specind in specsample:
    specdirectory = 'Spectra/'+specnames[specind]
    #print(specdirectory)

    fitdata = fits.getdata(specdirectory,ext=2)#import fits image

    metasize = len(fitdata[0])
    #print(metasize)
    if metasize == 126:
        redshift[specind] = fitdata[0][63]
        snmedian[specind] = np.max(fitdata[0][84])
        #for some reason they take the value from the 3rd filter
    else:
        redshift[specind] = fitdata[0][38]
        snmedian[specind] = np.median(fitdata[0][58])
    fitdata=0
print(snmedian)
