import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits

plt.style.use('mystyle')
plt.rcParams['lines.markersize'] = 2

posdata = fits.getdata('PositionsTable.fits',ext=1)
poslen = len(posdata)
raQSO = np.zeros(poslen)
decQSO = np.zeros(poslen)

for i in range(0,poslen):
    raQSO[i] = posdata[i][1]
    decQSO[i] = posdata[i][2]
# print(raQSO)
raQSO = (raQSO*np.pi)/180
decQSO = (decQSO*np.pi)/180

CARLAdata = fits.getdata('CARLA/CARLA_table_Aug_2013.fits',ext=1)
CARLAlen = len(CARLAdata)
raGC = np.zeros(CARLAlen)
decGC = np.zeros(CARLAlen)

for i in range(0,CARLAlen):
    raGC[i] = CARLAdata[i][1]
    decGC[i] = CARLAdata[i][2]
raGC = (raGC*np.pi)/180
decGC = (decGC*np.pi)/180
# r = np.sqrt(raGC**2+decGC**2)

plt.figure('CARLA Targets')
plt.subplot(111, projection="aitoff")
plt.plot(raGC,decGC,'r.')
plt.grid()

plt.figure('SDSS Targets')
plt.subplot(111, projection="aitoff")
plt.plot(raQSO,decQSO,'g.')
plt.grid()

plt.figure('Combined Targets')
plt.subplot(111, projection="aitoff")
plt.plot(raGC,decGC,'ro',markersize=6,mfc='none')
plt.plot(raGC,decGC,'k*')
plt.plot(raQSO,decQSO,'g.')
plt.grid()
# plt.colorbar()

plt.show()
