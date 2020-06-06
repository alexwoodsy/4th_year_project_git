import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate
import os

antidata = fits.getdata('Anti-Match-4point5degrees-z2point1.fits',ext=1)
antilen = len(antidata)
plate = np.zeros(antilen)
mjd = np.zeros(antilen)
fiberid = np.zeros(antilen)
print(antilen)

for i in range(0,antilen):
    plate[i] = antidata[i][4].astype(int)
    mjd[i] = antidata[i][5].astype(int)
    fiberid[i] = antidata[i][6].astype(int)
# print(plate)
#
f = open("download_output.txt", "a")

for j in range(0,antilen):
    string = '/spec-'
    print(str(int(plate[j]))+string+str(int(plate[j]))+'-'+str(int(mjd[j]))+'-'+str(int(fiberid[j])).zfill(4)+'.fits',file=f)

     # "-" mjd[j] "-" fiberid[j] ".fits")
f.close()
