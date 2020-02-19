import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits

<<<<<<< HEAD
data = fits.getdata('downloadspec/downloads/spec-7820-56984-0276.fits',ext=1)#import fits image
data2 = fits.getdata('downloadspec/downloads/spec-7820-56984-0276.fits',ext=2)#import fits image
print(data2)

# data2 = np.asarray(data2)
# print(data2)
# zindex = data2.find(2.4729445)
zindex = 38
z = data2[0][zindex]
print(len(data2[0]))

pred = (z+1)*1215.67
print(pred)
wlim = 1500
=======
data = fits.getdata('Spectra/spec-0834-52316-0243.fits',ext=1)#import fits image
wlim = 3000
>>>>>>> 10c55d4dbe0e141b11530a63bcb74c7d460ace7b
flux = np.zeros(len(data))
wlen = np.zeros(len(data))
model = np.zeros(len(data))
ivar = np.zeros(len(data))

for i in range(0,len(data)):
 flux[i] = data[i][0]
 ivar[i] = data[i][2]
 model[i] = data[i][7]
 wlen[i] = 10**(data[i][1])

<<<<<<< HEAD
out = flux#-model
=======
std = (ivar)**-0.5
ston = np.median(flux/std)
print("S/N ratio = ",ston)
>>>>>>> 10c55d4dbe0e141b11530a63bcb74c7d460ace7b

plt.plot(wlen[0:wlim],flux[0:wlim])
#plt.plot(wlen[0:wlim],model[0:wlim])
plt.xlabel('wavelength (Angstroms)')
plt.ylabel('Flux')

plt.show()
