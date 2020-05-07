import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate
import os

# antidata = fits.getdata('Anti-Match-4point5degrees-z2point1.fits',ext=1)
# antilen = len(antidata)

specnames = next(os.walk("C:/Users/jason/OneDrive - The University of Nottingham/4th Year/Project Files/Sky Spectra (Unused)"))[2]
spectot = len(specnames)

print(specnames)
plate = np.zeros(spectot)
mjd = np.zeros(spectot)
fiberid = np.zeros(spectot)

for i in range(0,spectot):
    plate[i] = specnames[i][5:9]
    mjd[i] = specnames[i][10:15]
    fiberid[i] = specnames[i][16:20]

# print(plate)
# print(mjd)
# print(fiberid)

f = open("metafix_output.txt", "a")
#
for j in range(0,spectot):
    print(str(int(plate[j])).zfill(4)+', '+str(int(mjd[j])).zfill(5)+', '+str(int(fiberid[j])).zfill(4),file=f)

f.close()
