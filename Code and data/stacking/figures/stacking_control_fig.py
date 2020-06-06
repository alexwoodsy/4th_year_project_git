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

#read in stacked data:
runsavename = 'lees_ALLbin'

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


#read in controll stack data:
inname = 'stacking/figures/Stacking data/' + runsavename+ '(control).fits'

controlstackdata = fits.getdata(inname,ext=1)
controlqsodata =fits.getdata(inname,ext=2)
controlstacklen = len(controlstackdata)
controlqsonumber = len(controlqsodata)


controlwlen = controlstackdata.field(0)
controlvrel = controlstackdata.field(1)
controlmed = controlstackdata.field(2)
controlmean = controlstackdata.field(3)

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
line0c, = ax0.plot(controlwlen, controlmean, color='g')
for i in range(1,12): #plot emission lines for reference
    ax0.plot(wlenline[i],mean[findval(wlen,wlenline[i])],'.')
ax0.set_ylim([-0.5, 1.5])
ax0.set_xlim([1100, 1300])
#the second subplot
ax1 = plt.subplot(gs[1], sharex = ax0, sharey = ax0)
line1, = ax1.plot(wlen, med, color='b')
line1c, = ax1.plot(controlwlen, controlmed, color='g')
for i in range(1,12): #plot emission lines for reference
    ax1.plot(wlenline[i],med[findval(wlen,wlenline[i])],'.',label = linename[i])
plt.setp(ax0.get_xticklabels(), visible=False)
# remove last tick label for the second subplot
yticks = ax1.yaxis.get_major_ticks()
yticks[-1].label1.set_visible(False)
# put leg on first subplot
ax0.legend((line0, line1, line0c), ('mean', 'median','control'), loc='lower left')
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
line0c, = ax0.plot(controlvrel, controlmean, color='g')
ax0.set_ylim([0, 1.5])
ax0.set_xlim([-3000, 3000])

#the second subplot
ax1 = plt.subplot(gs[1], sharex = ax0, sharey = ax0)
line1, = ax1.plot(vrel, med, color='b')
line1c, = ax1.plot(controlvrel, controlmed, color='g')
plt.setp(ax0.get_xticklabels(), visible=False)
# remove last tick label for the second subplot
yticks = ax1.yaxis.get_major_ticks()
yticks[-1].label1.set_visible(False)
# put leg on first subplot
ax0.legend((line0, line1, line0c), ('mean', 'median','control'), loc='lower left')

#labels
# Set common labels
fig.text(0.5, 0.04, r'$\delta$v ($kms^{-1})$', ha='center', va='center')
fig.text(0.06, 0.5, r'$F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$', ha='center', va='center', rotation='vertical')

# remove vertical gap between subplots
plt.subplots_adjust(hspace=.0)


#vrel binned pixels:
binsize = 100
vrelbins = np.arange(np.around(vrel[0],-1), np.around(vrel[-1],-1),binsize)
vrelbins = vrelbins[1:]
binind = np.zeros(len(vrelbins)).astype(int)
meanbinned = np.zeros(len(vrelbins))
medbinned = np.zeros(len(vrelbins))
controlmeanbinned = np.zeros(len(vrelbins))
controlmedbinned = np.zeros(len(vrelbins))

step =0
for i in range(0,len(vrelbins)):
    binind = findval(vrel, vrelbins[i])
    meanbinned[i] = np.mean(mean[step:binind])
    medbinned[i] = np.median(med[step:binind])
    controlmeanbinned[i] = np.mean(controlmean[step:binind])
    controlmedbinned[i] = np.median(controlmed[step:binind])
    step = binind

fig = plt.figure('vrel_binned')
# set height ratios for sublots
gs = plt.GridSpec(2, 1, height_ratios=[1, 1])

# the fisrt subplot
ax0 = plt.subplot(gs[0])
#ax0.set_yscale("log")
line0, = ax0.step(vrelbins, meanbinned, color='r')
line0c, = ax0.step(vrelbins, controlmeanbinned, color='g')
ax0.set_ylim([0, 1.5])
ax0.set_xlim([-3000, 3000])

#the second subplot
ax1 = plt.subplot(gs[1], sharex = ax0, sharey = ax0)
line1, = ax1.step(vrelbins, medbinned, color='b')
line1c, = ax1.step(vrelbins, controlmedbinned, color='g')
plt.setp(ax0.get_xticklabels(), visible=False)
# remove last tick label for the second subplot
yticks = ax1.yaxis.get_major_ticks()
yticks[-1].label1.set_visible(False)
# put leg on first subplot
ax0.legend((line0, line1, line0c), ('mean', 'median','control'), loc='lower left')

#labels
# Set common labels
fig.text(0.5, 0.04, r'$\delta$v ($kms^{-1})$', ha='center', va='center')
fig.text(0.06, 0.5, r'$F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$', ha='center', va='center', rotation='vertical')

# remove vertical gap between subplots
plt.subplots_adjust(hspace=.0)

###absorbtion line fitting####
#fit controll level
def polynomial(x, a, b):
    return a*(x)+b

contdatarange = np.arange(findval(vrelbins,-1000), findval(vrelbins,1000))
contpopt, pcov = cf(polynomial, vrelbins[contdatarange], controlmeanbinned[contdatarange], bounds =([-100,-1000],[100,1000]))
conta, contb = contpopt
control_fit = polynomial(vrelbins[contdatarange], *contpopt)
control_level = 1


#fitting needs initial data so extract data about abs line
def guassian(x, amp, mean, std):
    return amp*np.exp(-((x-mean)**2)/(2*(std)**2))

figabs, ax = plt.subplots(1,1,num='Absorption line fitting')
datarange = np.arange(findval(vrelbins,-1000), findval(vrelbins,1000))
popt, pcov = cf(guassian, vrelbins[datarange], control_fit-meanbinned[datarange], bounds =([-1,-700,0],[1.5,700,1000]))
meanamp, meanmean, meanstd = popt
ax.step(vrelbins, meanbinned,label='stack data')
ax.plot(vrelbins[datarange], control_fit-guassian(vrelbins[datarange], *popt), 'r-',label='fitting parmaters: amp=%5.3f, mean=%5.3f, std=%5.3f' % tuple(popt))

ax.step(vrelbins, controlmeanbinned,label='control data')
ax.plot(vrelbins[contdatarange], control_fit, 'b-',label='fitting parmaters: a=%5.3f, b=%5.3f' % tuple(contpopt))

ax.set_xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
ax.set_ylabel(r'$<F>$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')
ax.legend()


plt.show()
