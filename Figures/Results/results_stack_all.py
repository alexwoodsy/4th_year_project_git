import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
from scipy.optimize import curve_fit as cf
import os, time
#import fitting
import sys

plt.style.use('mystyle')
def findval(array,val):
    array = np.asarray(array)
    ind = np.abs(array - val).argmin()
    return ind

runs = ['ours_ALLbin','old stacking data/lee_Allcarla']

#read in stacked data:
for runsavename in runs:
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
    print(controlqsonumber)

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

    if runsavename == runs[0]:
        oursvrelbins = vrelbins
        oursmeanbinned = meanbinned
        oursmedbinned = medbinned
        ourscontrolmeanbinned = controlmeanbinned
        ourscontrolmedbinned = controlmedbinned
    elif runsavename == runs[1]:
        leesvrelbins = vrelbins
        leesmeanbinned = meanbinned
        leesmedbinned = medbinned
        leescontrolmeanbinned = controlmeanbinned
        leescontrolmedbinned = controlmedbinned



fig = plt.figure('vrel_binned')
# set height ratios for sublots
gs = plt.GridSpec(2, 2, height_ratios=[1, 1])

# the fisrt subplot
ax0 = plt.subplot(gs[0])
#ax0.set_yscale("log")
ax0.step(oursvrelbins, oursmeanbinned, color='#ff0000')
ax0.step(oursvrelbins, ourscontrolmeanbinned, color='#00ff00')
ax0.set_ylim([0, 1.5])
ax0.set_xlim([-3000, 3000])

#the second subplot
ax1 = plt.subplot(gs[2], sharex = ax0, sharey = ax0)
ax1.step(oursvrelbins, oursmedbinned, color='#0000ff')
ax1.step(oursvrelbins, ourscontrolmedbinned, color='#00ff00')
ax1.set_xlabel(r'$\delta$v ($kms^{-1})$')
plt.setp(ax0.get_xticklabels(), visible=False)
# remove last tick label for the second subplot
yticks = ax1.yaxis.get_major_ticks()
yticks[-1].label1.set_visible(False)

# the 3rd subplot
ax2 = plt.subplot(gs[1], sharex = ax0)
#ax0.set_yscale("log")
ax2.step(leesvrelbins, leesmeanbinned, color='#ff0000')
ax2.step(leesvrelbins, leescontrolmeanbinned, color='#00ff00')
ax2.set_ylim([0, 1.5])
ax2.set_xlim([-3000, 3000])

#the 4th subplot
ax3 = plt.subplot(gs[3], sharex = ax0, sharey = ax2)
ax3.step(leesvrelbins, leesmedbinned, color='#0000ff')
ax3.step(leesvrelbins, leescontrolmedbinned, color='#00ff00')
ax3.set_xlabel(r'$\delta$v ($kms^{-1})$')
plt.setp(ax2.get_xticklabels(), visible=False)
# remove last tick label for the second subplot
yticks = ax3.yaxis.get_major_ticks()
yticks[-1].label1.set_visible(False)

fig.text(0.05, 0.5, r'$\tilde{F}$', ha='center', va='center', rotation='vertical')

# remove vertical gap between subplots
plt.subplots_adjust(hspace=.0)
plt.show()
