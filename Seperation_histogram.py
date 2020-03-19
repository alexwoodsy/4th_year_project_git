import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
import os

# plt.style.use('mystyle')

positions = fits.getdata('HalfDegreeSeperation.fits',ext=1)
poslen = len(positions)
print(poslen)
seperationsq = np.zeros(poslen)

for j in range(0):
    seperationsq[j] = positions[j][116]

print(seperationsq)
#
plt.figure()
# logsep = np.log(seperation*seperation)
logsep = seperationsq
step = 0.5
binning = np.arange(min(logsep)-1, max(logsep) + 1, step)
n, bins, patches = plt.hist(x = logsep, bins = binning, alpha=1, rwidth=1)
print(n)
#plt.grid(axis='y', alpha=0.4)
plt.xlabel(r'$Seperation^2$')
# plt.xticks(binning)
plt.ylabel('N')

plt.show()
