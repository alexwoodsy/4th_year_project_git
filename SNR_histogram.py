import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate
import os

specnames = next(os.walk('Spectra'))[2] #dir is your directory path as string
spectot = len(specnames)

#add indexing for epctra in file to allow loop over all
#specind=1
specsample = np.array([spectot-1])
# print(specsample)
totivar = []
totflux = []

for specind in range(0,max(specsample)):
    specdirectory = 'Spectra/'+specnames[specind]
    # print(specdirectory)
    data = fits.getdata(specdirectory,ext=1)#import fits image
    # speclen = len(data)
    # flux = np.zeros(speclen)
    # # wlen = np.zeros(speclen)
    # # model = np.zeros(speclen)
    # ivar = np.zeros(speclen)

    medflux = np.median(data[0])
    medivar = np.median(data[2])

    # for i in range(0,speclen):
    #      flux[i] = data[i][0]
    #      ivar[i] = data[i][2]
    #      # model[i] = data[i][7]
    #      # wlen[i] = 10**(data[i][1])
    # medivar = np.median(ivar)
    # medflux = np.median(flux)
    # print(medivar)
    totivar.append(medivar)
    totflux.append(medflux)

print(totivar)
#ston calculation
std = (np.asarray(totivar))**-0.5
ston = np.asarray(totflux)/std
# print(ston)
# print(len(ivar))
logsn = np.log(ston)
binning = np.arange(min(logsn), max(logsn) + 1, 1)

n, bins, patches = plt.hist(x = logsn, bins = binning,color='green', alpha=0.7, rwidth=0.9)
# print(n)
plt.grid(axis='y', alpha=0.4)
plt.xlabel('SNR')
plt.ylabel('Frequency')
plt.title('SNR Histogram')

plt.show()
