import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits

data = fits.getdata('Spectra/spec-0834-52316-0243.fits',ext=1)#import fits image
wlim = 3000
flux = np.zeros(len(data))
wlen = np.zeros(len(data))
model = np.zeros(len(data))
ivar = np.zeros(len(data))

for i in range(0,len(data)):
 flux[i] = data[i][0]
 ivar[i] = data[i][2]
 model[i] = data[i][7]
 wlen[i] = 10**(data[i][1])

std = (ivar)**-0.5
ston = np.median(flux/std)
print("S/N ratio = ",ston)

plt.plot(wlen[0:wlim],flux[0:wlim])
#plt.plot(wlen[0:wlim],model[0:wlim])
plt.xlabel('wavelength (Angstroms)')
plt.ylabel('Flux')

plt.show()
