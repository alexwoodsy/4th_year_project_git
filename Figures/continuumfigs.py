import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
from scipy.optimize import curve_fit as cf
import os, time
#import fitting
import sys
plt.style.use('mystyle')
# Optionally set font to Computer Modern to avoid common missing font errors


#mpl.rcParams['text.latex.preamble'] = [r'\renewcommand{\seriesdefault}{\bfdefault}']
#mpl.rcParams['text.latex.preamble'] = [r'\renewcommand{\familydefault}{cmr}']
#mpl.rcParams['text.latex.preamble'] = [r'\boldmath']
#mpl.rcParams['text.latex.preamble'] = [r'\symrm']


#imports the  FITTED spectra from the spectra folder
specnames = next(os.walk('Fitted Spectra'))[2]
spectot = len(specnames)

#import old to get sdss name:
oldspecnames = next(os.walk('Spectra'))[2]
oldspectot = len(specnames)




speccut =False


#import spectrum data fitted by us:
spec = specnames[0]

#get spec name from old data:
sdssdata = fits.getdata('PositionsTable.fits',ext=1)
i = 0
namecheck = False
for i in range(len(sdssdata)):
    file = 'spec-'+str(sdssdata[i][4]).zfill(4) +'-'+ str(sdssdata[i][5]).zfill(5)+'-' + str(sdssdata[i][6]).zfill(4)
    if str(file) == spec[0:20]:
        sdssspecname = sdssdata[i][0]



specdirectory = 'Fitted Spectra/'+spec
data = fits.getdata(specdirectory,ext=1) #read in ours
speclen = len(data)
#get spec/fit data
wlen = data.field(0)
flux = data.field(1)
medcontfit = data.field(2)
maxcontfit = data.field(3)
medcontfit_2 = data.field(4)
maxcontfit_2 = data.field(5)
max_medcontfit = data.field(6)
contfit = max_medcontfit

metadata = fits.getdata(specdirectory,ext=2)
z = np.round(metadata[0][0],3)


if speccut == True: #cut spec to same length as lee (1030 to 1600 in rest frame)
    cutstart = findval(wlen,(1030*(1+redshift)))
    cutend = findval(wlen,(1600*(1+redshift)))
    if cutend != 0 and cutstart != 0:
        cut = np.arange(cutstart,cutend)
        contfit = contfit[cut]
        flux = flux[cut]
        wlen = wlen[cut]


#continuum normalisation
normspec = flux/contfit

fig = plt.figure('cont')
# set height ratios for sublots
gs = plt.GridSpec(3, 1, height_ratios=[1, 1, 1])
# the fisrt subplot
ax0 = plt.subplot(gs[0])
#ax0.set_yscale("log")
line0, = ax0.plot(wlen, flux)
name0 = ax0.text(8750,140,sdssspecname)
z0 = ax0.text(9500,120,r'z = '+str(z))

#second
spec = specnames[10]

#get spec name from old data:
sdssdata = fits.getdata('PositionsTable.fits',ext=1)
i = 0
namecheck = False
for i in range(len(sdssdata)):
    file = 'spec-'+str(sdssdata[i][4]).zfill(4) +'-'+ str(sdssdata[i][5]).zfill(5)+'-' + str(sdssdata[i][6]).zfill(4)
    if str(file) == spec[0:20]:
        sdssspecname = sdssdata[i][0]



specdirectory = 'Fitted Spectra/'+spec
data = fits.getdata(specdirectory,ext=1) #read in ours
speclen = len(data)
#get spec/fit data
wlen = data.field(0)
flux = data.field(1)
medcontfit = data.field(2)
maxcontfit = data.field(3)
medcontfit_2 = data.field(4)
maxcontfit_2 = data.field(5)
max_medcontfit = data.field(6)
contfit = max_medcontfit

metadata = fits.getdata(specdirectory,ext=2)
z = np.round(metadata[0][0],3)


if speccut == True: #cut spec to same length as lee (1030 to 1600 in rest frame)
    cutstart = findval(wlen,(1030*(1+redshift)))
    cutend = findval(wlen,(1600*(1+redshift)))
    if cutend != 0 and cutstart != 0:
        cut = np.arange(cutstart,cutend)
        contfit = contfit[cut]
        flux = flux[cut]
        wlen = wlen[cut]


#continuum normalisation
normspec = flux/contfit

# the fisrt subplot
ax1 = plt.subplot(gs[1], sharex = ax0)
#ax0.set_yscale("log")
line1, = ax1.plot(wlen, flux)
name1 = ax1.text(8665 ,11,sdssspecname)
z1 = ax1.text(9500,8.5,r'z = '+str(z))




#the second subplot

spec = specnames[1000]

#get spec name from old data:
sdssdata = fits.getdata('PositionsTable.fits',ext=1)
i = 0
namecheck = False
for i in range(len(sdssdata)):
    file = 'spec-'+str(sdssdata[i][4]).zfill(4) +'-'+ str(sdssdata[i][5]).zfill(5)+'-' + str(sdssdata[i][6]).zfill(4)
    if str(file) == spec[0:20]:
        sdssspecname = sdssdata[i][0]



specdirectory = 'Fitted Spectra/'+spec
data = fits.getdata(specdirectory,ext=1) #read in ours
speclen = len(data)
#get spec/fit data
wlen = data.field(0)
flux = data.field(1)
medcontfit = data.field(2)
maxcontfit = data.field(3)
medcontfit_2 = data.field(4)
maxcontfit_2 = data.field(5)
max_medcontfit = data.field(6)
contfit = max_medcontfit

metadata = fits.getdata(specdirectory,ext=3)
z = np.round(metadata[0][0],3)


if speccut == True: #cut spec to same length as lee (1030 to 1600 in rest frame)
    cutstart = findval(wlen,(1030*(1+redshift)))
    cutend = findval(wlen,(1600*(1+redshift)))
    if cutend != 0 and cutstart != 0:
        cut = np.arange(cutstart,cutend)
        contfit = contfit[cut]
        flux = flux[cut]
        wlen = wlen[cut]


#continuum normalisation
normspec = flux/contfit


ax2 = plt.subplot(gs[2], sharex = ax0)
line2, = ax2.plot(wlen, flux)
name2 = ax2.text(8665,25,sdssspecname)
z2 = ax2.text(9500,21,r'z = '+str(z))
# remove last tick label for the second subplot
yticks = ax1.yaxis.get_major_ticks()
yticks[-1].label1.set_visible(False)
# put leg on first subplot
#ax0.legend((line0, line1), ('mean', 'median'), loc='lower left')

#labels
# Set common labels
fig.text(0.5, 0.03, r'{$\lambda$ ($\mathrm{\AA}$)', ha='center', va='center')
fig.text(0.05, 0.5, r'$F$ $(10^{-17} $erg$ s^{-1}cm^{-2}\mathrm{\AA}^{-1})$', ha='center', va='center', rotation='vertical')
# remove vertical gap between subplots
plt.xlim(wlen[0], wlen[-1])
plt.subplots_adjust(hspace=.0)
plt.show()















#
