import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits

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
flux = np.zeros(len(data))
wlen = np.zeros(len(data))
model = np.zeros(len(data))

for i in range(0,len(data)):
 flux[i] = data[i][0]
 model[i] = data[i][7]
 wlen[i] = 10**(data[i][1])

out = flux#-model

plt.plot(wlen[0:wlim],out[0:wlim])
plt.xlabel('wavelength (Angstroms)')
plt.ylabel('Flux')

plt.show()
