import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy.interpolate import interp1d
import os

specnames = next(os.walk('Spectra'))[2] #dir is your directory path as string
spectot = len(specnames)

#add indexing for epctra in file to allow loop over all
i=1
specdirectory = 'Spectra/'+specnames[i]
print(specdirectory)

data = fits.getdata(specdirectory,ext=1)#import fits image
wlim = len(data)
flux = np.zeros(len(data))
wlen = np.zeros(len(data))
model = np.zeros(len(data))
ivar = np.zeros(len(data))

for i in range(0,len(data)):
 flux[i] = data[i][0]
 ivar[i] = data[i][2]
 model[i] = data[i][7]
 wlen[i] = 10**(data[i][1])

#ston calculation
std = (ivar)**-0.5
ston = np.median(flux/std)
print("S/N ratio = ",ston)


#cont fitting
intervals = 30
window = int(len(data)/intervals)
step = 0
i = 0
intervalwlen = np.zeros(intervals)
winpeak = np.zeros(intervals)

while (step+window) <= len(data):
    windata = flux[step:(step+window)]
    winpeakind = step + np.argmax(windata)
    winpeak[i] = np.max(windata)
    intervalwlen[i] = wlen[winpeakind]
    step = step + window
    i = i + 1

intpol = interp1d(intervalwlen, winpeak, kind='cubic')
xnew = np.linspace(intervalwlen[0],intervalwlen[-1], num=100, endpoint=True)


#plotting

plt.plot(wlen[0:wlim],flux[0:wlim])
plt.plot(intervalwlen[0:wlim],winpeak[0:wlim],'*')
plt.plot(xnew,intpol(xnew),'--')
plt.xlabel('wavelength (Angstroms)')
plt.ylabel('Flux')

plt.show()
