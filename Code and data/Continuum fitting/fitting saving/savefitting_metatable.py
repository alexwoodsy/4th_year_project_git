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

specfilename = []
redshift = []
stonall = []
stonforest = []
lyalpha = []
lyalphaind = []
gcname = []
gcredshift = []
gc_qso_sep = []
gclyalpha = []
gclyalphaind = []
stackmsg = []
leecheckcol = []
leecontflag = []


maxlenspecname = 0
for spec in specnames:
    hdul = fits.open('Fitted Spectra/'+spec)
    hdunum = len(hdul)
    if hdunum == 3: #reads correct thing in if lee data present
        fitdata = fits.getdata('Fitted Spectra/'+spec, ext=2)
    else:
        fitdata = fits.getdata('Fitted Spectra/'+spec, ext=3)

    specfilename = np.append(specfilename, spec)
    if len(spec) > maxlenspecname:
        maxlenspecname = len(spec)
    redshift = np.append(redshift, fitdata[0][0])
    stonall = np.append(stonall, fitdata[0][1])
    stonforest = np.append(stonforest, fitdata[0][2])
    lyalpha = np.append(lyalpha, fitdata[0][3])
    lyalphaind = np.append(lyalphaind, fitdata[0][4])
    gcname = np.append(gcname, fitdata[0][5])
    gcredshift = np.append(gcredshift, fitdata[0][6])
    gc_qso_sep = np.append(gc_qso_sep, fitdata[0][7])
    gclyalpha = np.append(gclyalpha, fitdata[0][8])
    gclyalphaind = np.append(gclyalphaind, fitdata[0][9])
    stackmsg = np.append(stackmsg, fitdata[0][10])
    leecheckcol = np.append(leecheckcol, fitdata[0][11])
    leecontflag = np.append(leecontflag, fitdata[0][12])


#put in from_columns
specfilenamecol = fits.Column(name='SPEC_FILE_NAME', array = specfilename, format=str(maxlenspecname)+'A')
redshiftcol = fits.Column(name='QSO_Z', array = redshift, format='F')
stonallcol = fits.Column(name='STON_ALL', array = stonall, format='F')
stonforestcol = fits.Column(name='STON_FOREST', array= stonforest, format='F')
lyalphacol = fits.Column(name='QSO_LYa', array = lyalpha, format='F')
lyalphaindcol = fits.Column(name='QSO_LYa_INDEX', array = lyalphaind, format='K')
gcnamecol = fits.Column(name='GC_NAME', array = gcname, format='20A')
gcredshiftcol = fits.Column(name='GC_Z', array = gcredshift, format='F')
gc_qso_sepcol = fits.Column(name='GC_SEPERATION', array = gc_qso_sep, format='F')
gclyalphacol = fits.Column(name='GC_LYa', array = gclyalpha, format='F')
gclyalphaindcol = fits.Column(name='GC_LYa_INDEX', array = gclyalphaind, format='K')
stackmsgcol = fits.Column(name='STACK_STATUS', array = stackmsg, format='20A')
leecheckcol = fits.Column(name='LEE_CHECK', array = leecheckcol, format='K')
leecontflagcol = fits.Column(name='LEE_CONT_FLAG', array = leecontflag, format='K')

metadata = fits.BinTableHDU.from_columns([specfilenamecol, redshiftcol, stonallcol, stonforestcol, lyalphacol,
lyalphaindcol, gcnamecol, gcredshiftcol, gc_qso_sepcol, gclyalphacol, gclyalphaindcol, stackmsgcol,leecheckcol, leecontflagcol])

primary = fits.PrimaryHDU()
hdul = fits.HDUList([primary, metadata])

hdul.writeto('metatable_V2.fits',overwrite = True)

#################specfilenamecol,



# if len(spec) != 35: #standard length == 35 hence anything more == spec in 2 clusters
#     print(spec[21:-16])
#



#
