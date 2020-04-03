import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
import os
#import fitting
import sys
#sys.path.append('C:/Users/jason/GIT/4th_year_project_git/Continuum Fitting')
sys.path.append('C:/Users/alexw/Documents/GitHub/4th_year_project_git/Continuum fitting')


#imports the spectra from the spectra folder
specnames = next(os.walk('Fitted Spectra'))[2]
spectot = len(specnames)

nullerror = []
count = 0
for spec in specnames:
    specdirectory = 'Fitted Spectra/'+spec
    fitdata = fits.getdata(specdirectory,ext=2)#read in fits meta data
    redshift= fitdata[0][0]
    stonall = fitdata[0][1]
    stonforest = fitdata[0][2]
    lyalpha = fitdata[0][3]
    lyalphaind = fitdata[0][4]
    gcname = fitdata[0][5]
    gcredshift = fitdata[0][6]
    gc_qso_sep = fitdata[0][7]
    gclyalpha = fitdata[0][8]
    gclyalphaind = fitdata[0][9]
    stackmsg = fitdata[0][10]

    if gcname == 'null':
        count = count + 1
        print('found ' + str(count))
        nullerror.append(spec)
print(nullerror)
