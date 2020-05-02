#stacking v2
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
from scipy.optimize import curve_fit as cf
import os, time
#import fitting
import sys

#plt.style.use('mystyle') #path C:\Users\alexw\AppData\Local\Programs\Python\Python37\Lib\site-packages\matplotlib\mpl-data\stylelib


arr2d = np.zeros([9,9])
arr2d[0:9, 0:9] = np.nan

for i in range(0,8):
    arr2d[i,5:9] = np.arange(5,9)

arr2d[3, 1:5] = np.arange(1,5)
arr2d[4, 3:5] = np.arange(3,5)
print(arr2d)

med = np.nanmean(arr2d, axis=0)
print(med)
med = med[~np.isnan(med)]

print(med)
