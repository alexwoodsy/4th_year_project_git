import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
import os

plt.style.use('mystyle')

###############################################################################################
# issue with plt.hist causing massive runtime, will fix later but code is correct for this plot
###############################################################################################

positions = fits.getdata('HalfDegreeSeperation.fits',ext=1)
poslen = len(positions)
print(poslen)
# print(positions[0][116])
seperationsq = np.zeros(10)

for j in range(0,10):
    seperationsq[j] = positions[j][116]

print(seperationsq)
#
plt.figure()
# logsep = np.log(seperation*seperation)
logsep = (seperationsq)/3600
step = 0.5
binning = np.arange(min(logsep)-1, max(logsep) + 1, step)
print(step)
plt.hist(x = logsep, bins = binning, alpha=1, rwidth=1, histtype='step')
# print(n)
#plt.grid(axis='y', alpha=0.4)
plt.xlabel(r'$Seperation^2 (arcmin^2)$')
# plt.xticks(binning)
plt.ylabel('N')

plt.show()
