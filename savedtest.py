import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
import os
#import fitting
import sys
sys.path.append('C:/Users/jason/GIT/4th_year_project_git/Continuum Fitting')

# specdirectory = 'Spectra/spec-0411-51817-0290.fits'
#
# data = fits.getdata(specdirectory,ext=2)
#
# print(data)

n = np.ones([1,3])


for i in range(0,4):
    n2 = np.zeros([i,3])
    # print(n)

    n = np.append(n,n2,axis=0)
    print(n)

m = np.median(n,axis=0)
print(m)
