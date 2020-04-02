import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
import os
#import fitting
import sys
sys.path.append('C:/Users/jason/GIT/4th_year_project_git/Continuum Fitting')

specdirectory = 'Fitted Spectra/spec-0411-51817-0290-prefitted.fits'

data = fits.getdata(specdirectory,ext=0)

print(data)
