import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits

data = fits.getdata('Testing/spec-0266-51602-0050.fits',ext=1)#import fits image
print(data)
wlim = 1500
flux = np.zeros(len(data))
wlen = np.zeros(len(data))
model = np.zeros(len(data))

for i in range(0,len(data)):
 flux[i] = data[i][0]
 model[i] = data[i][7]
 wlen[i] = 10**(data[i][1])

out = flux-model

plt.plot(wlen[0:wlim],out[0:wlim])
plt.xlabel('wavelength (Angstroms)')
plt.ylabel('Flux')

plt.show()
