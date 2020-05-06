import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
from scipy.optimize import curve_fit as cf
import os, time
#import fitting
import sys

#plt.style.use('mystyle')
def findval(array,val):
    array = np.asarray(array)
    ind = np.abs(array - val).argmin()
    return ind

#read ion stacked data:
runsavename = 'max_med_Allcarla'
print(np.random.random(1))

inname = 'stacking/figures/Stacking data/' + runsavename+ '.fits'

stackdata = fits.getdata(inname,ext=1)
carladata = fits.getdata(inname,ext=2)
qsodata =fits.getdata(inname,ext=3)
stacklen = len(stackdata)
qsonumber = len(qsodata)
print(qsonumber)
carlanumber = len(carladata)

wlen = stackdata.field(0)
vrel = stackdata.field(1)
med = stackdata.field(2)
mean = stackdata.field(3)

#read in abs lines
refline = open('refline.txt', 'r')
linename = []
wlenline = np.array([])
for line in refline:
    line = line.strip()
    columns = line.split()
    wlenline = np.append(wlenline, float(columns[0]))
    linename.append(columns[1])
refline.close()
linename = linename[0:12] #restrict higher WLEN lines
wlenline = wlenline[0:12]



fig = plt.figure('mean and median')
# set height ratios for sublots
gs = plt.GridSpec(2, 1, height_ratios=[1, 1])
# the fisrt subplot
ax0 = plt.subplot(gs[0])
#ax0.set_yscale("log")
line0, = ax0.plot(wlen, mean, color='r')
for i in range(1,12): #plot emission lines for reference
    ax0.plot(wlenline[i],mean[findval(wlen,wlenline[i])],'.')
ax0.set_ylim([-0.5, 1.5])
ax0.set_xlim([1100, 1300])
#the second subplot
ax1 = plt.subplot(gs[1], sharex = ax0, sharey = ax0)
line1, = ax1.plot(wlen, med, color='b')
for i in range(1,12): #plot emission lines for reference
    ax1.plot(wlenline[i],med[findval(wlen,wlenline[i])],'.',label = linename[i])
plt.setp(ax0.get_xticklabels(), visible=False)
# remove last tick label for the second subplot
yticks = ax1.yaxis.get_major_ticks()
yticks[-1].label1.set_visible(False)
# put leg on first subplot
ax0.legend((line0, line1), ('mean', 'median'), loc='lower left')
ax1.legend(loc='lower right')
#labels
# Set common labels
fig.text(0.5, 0.04, r'$\lambda$ ($\mathrm{\AA}$)', ha='center', va='center')
fig.text(0.06, 0.5, r'$F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$', ha='center', va='center', rotation='vertical')
# remove vertical gap between subplots
plt.subplots_adjust(hspace=.0)

fig = plt.figure('vrel')
# set height ratios for sublots
gs = plt.GridSpec(2, 1, height_ratios=[1, 1])

# the fisrt subplot
ax0 = plt.subplot(gs[0])
#ax0.set_yscale("log")
line0, = ax0.plot(vrel, mean, color='r')
ax0.set_ylim([0, 1.5])
ax0.set_xlim([-3000, 3000])

#the second subplot
ax1 = plt.subplot(gs[1], sharex = ax0, sharey = ax0)
line1, = ax1.plot(vrel, med, color='b')
plt.setp(ax0.get_xticklabels(), visible=False)
# remove last tick label for the second subplot
yticks = ax1.yaxis.get_major_ticks()
yticks[-1].label1.set_visible(False)
# put leg on first subplot
ax0.legend((line0, line1), ('mean', 'median'), loc='lower left')

#labels
# Set common labels
fig.text(0.5, 0.04, r'$\delta$v ($kms^{-1})$', ha='center', va='center')
fig.text(0.06, 0.5, r'$F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$', ha='center', va='center', rotation='vertical')


# remove vertical gap between subplots
plt.subplots_adjust(hspace=.0)




plt.show()
