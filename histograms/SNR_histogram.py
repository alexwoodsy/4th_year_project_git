import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate
import os

plt.style.use('mystylehistogram')

specnames = next(os.walk('Spectra'))[2] #dir is your directory path as string
spectot = len(specnames)

#add indexing for epctra in file to allow loop over all
#specind=1
specsample = np.arange(0,spectot-1)
# specsample = np.arange(0,200)
# print(specsample)
totivar = []
totflux = []

for specind in specsample:
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

#print(totivar)
#ston calculation
std = (np.asarray(totivar))**-0.5
ston = np.asarray(totflux)/std
# print(ston)
# print(len(ivar))
logsn = np.log(ston)
step = 0.5
binning = np.arange(int(min(logsn)-step)-0.5, int(max(logsn)+step) + 1, step)
n, bins, patches = plt.hist(x = logsn, bins = binning, alpha=1, rwidth=1)
# print(n)
#plt.grid(axis='y', alpha=0.4)
plt.xlabel(r'$log(S/N_{Ly\alpha})$')
plt.xticks(binning)

plt.ylabel('N')
#plt.title('SNR Histogram')

plt.show()
