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

#plot initialisation:
fig = plt.figure('vrel_binned')
# set height ratios for sublots
gs = plt.GridSpec(2, 3, height_ratios=[1, 1])


plotcounter = 0
namecounter = 0

runs = ['ours_200to400bin','ours_200to400bin_sn_2','ours_200to400bin_sn_1.5']
panelname = [r'S/N $>$ 4',r'S/N $>$ 2',r'S/N $>$ 1.5']
#read in stacked data:
for runsavename in runs:
    inname = 'stacking/figures/Stacking data/' + runsavename+ '.fits'

    stackdata = fits.getdata(inname,ext=1)
    carladata = fits.getdata(inname,ext=2)
    qsodata =fits.getdata(inname,ext=3)
    stacklen = len(stackdata)
    qsonumber = len(qsodata)

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
    linename = linename[2:12] #restrict higher WLEN lines
    wlenline = wlenline[2:12]
    c = 299792.458
    lam = wlenline
    lam_em = 1215.67
    vrelline = c*((lam  - lam_em)/lam_em)





    #vrel binned pixels:
    binsize = 100
    vrelbins = np.arange(np.around(vrel[0],-1), np.around(vrel[-1],-1),binsize)
    vrelbins = vrelbins[1:]
    binind = np.zeros(len(vrelbins)).astype(int)
    meanbinned = np.zeros(len(vrelbins))
    medbinned = np.zeros(len(vrelbins))
    controlmeanbinned = np.zeros(len(vrelbins))
    controlmedbinned = np.zeros(len(vrelbins))

    print(runsavename)
    print('qso number  = '+str(qsonumber))
    print('carla number  = '+str(carlanumber))
    print('control number  = '+str(controlqsonumber))

    step =0
    for i in range(0,len(vrelbins)):
        binind = findval(vrel, vrelbins[i])
        meanbinned[i] = np.mean(mean[step:binind])
        medbinned[i] = np.median(med[step:binind])
        controlmeanbinned[i] = np.mean(controlmean[step:binind])
        controlmedbinned[i] = np.median(controlmed[step:binind])
        step = binind

    # the fisrt subplot seperately
    if plotcounter == 0:
        axl0 = plt.subplot(gs[plotcounter],label=runsavename)
        #for line in vrelline: #plot em lines
        #    axl0.plot(np.array([line,line]),np.array([-1,5]),'#00ffff', linestyle =':')
        #ax0.set_yscale("log")
        axl0.plot(np.array([vrelbins[0],vrelbins[-1]]),np.array([1,1]),'#FFA500', linestyle =':')
        axl0.plot(np.array([1,1]), np.array([-1,5]),'#FFA500', linestyle =':')
        axl0.step(vrelbins, meanbinned, color='#ff0000')
        axl0.step(vrelbins, controlmeanbinned, color='#00ff00')
        axl0.text(-3800, 0.52,panelname[namecounter])



        axr0 = plt.subplot(gs[plotcounter+3],sharex = axl0, sharey = axl0, label=runsavename)
        #ax0.set_yscale("log")
        #for line in vrelline: #plot em lines
        #    axr0.plot(np.array([line,line]),np.array([-1,5]),'#00ffff', linestyle =':')
        axr0.plot(np.array([vrelbins[0],vrelbins[-1]]),np.array([1,1]),'#FFA500', linestyle =':')
        axr0.plot(np.array([1,1]), np.array([-1,5]),'#FFA500', linestyle =':')
        axr0.step(vrelbins, medbinned, color='#0000ff')
        axr0.step(vrelbins, controlmedbinned, color='#00ff00')
        axr0.set_ylim([0.55, 1.05])
        axr0.set_xlim([-2500, 3000])
        axr0.text(-4000, 0.52,'$N_{QSO}^{LOS}$ = '+str(qsonumber))
        axr0.text(-4000, 0.6,'$N_{CARLA}$ = '+str(carlanumber))
        axr0.text(600, 0.52,'$N_{QSO}^{Control}$ = '+str(controlqsonumber))
        axr0.set_xlabel(r'$\delta v$ (km s$^{-1})$')
    else:
        axl = plt.subplot(gs[plotcounter],sharex = axl0, sharey = axl0, label=runsavename)
        #ax0.set_yscale("log")
        #for line in vrelline: #plot em lines
        #    axl.plot(np.array([line,line]),np.array([-1,5]),'#00ffff', linestyle =':')
        axl.plot(np.array([vrelbins[0],vrelbins[-1]]),np.array([1,1]),'#FFA500', linestyle =':')
        axl.plot(np.array([1,1]), np.array([-1,5]),'#FFA500', linestyle =':')
        axl.step(vrelbins, meanbinned, color='#ff0000')
        axl.step(vrelbins, controlmeanbinned, color='#00ff00')
        axl.text(-3800, 0.52,panelname[namecounter])
        if runsavename != runs[-1]: #
            yticks = axl.yaxis.get_major_ticks()
            xticks = axl.xaxis.get_major_ticks()
            yticks[-1].label1.set_visible(False)
            xticks[-1].label1.set_visible(False)

        axr = plt.subplot(gs[plotcounter+3],sharex = axl0, sharey = axl0, label=runsavename)
        #ax0.set_yscale("log")
        #for line in vrelline: #plot em lines
        #    axr.plot(np.array([line,line]),np.array([-1,5]),'#00ffff', linestyle =':')
        axr.plot(np.array([vrelbins[0],vrelbins[-1]]),np.array([1,1]),'#FFA500', linestyle =':')
        axr.plot(np.array([1,1]), np.array([-1,5]),'#FFA500', linestyle =':')
        axr.step(vrelbins, medbinned, color='#0000ff')
        axr.step(vrelbins, controlmedbinned, color='#00ff00')
        axr.text(-4000, 0.52,'$N_{QSO}^{LOS}$ = '+str(qsonumber))
        axr.text(-4000, 0.6,'$N_{CARLA}$ = '+str(carlanumber))
        axr.text(300, 0.52,'$N_{QSO}^{Control}$ = '+str(controlqsonumber))
        axr.set_ylim([0.55, 1.05])
        axr.set_xlim([-2500, 3000])
        if runsavename != runs[-1]: #
            yticks = axr.yaxis.get_major_ticks()
            xticks = axr.xaxis.get_major_ticks()
            yticks[-1].label1.set_visible(False)
            xticks[-1].label1.set_visible(False)


        axr.set_xlabel(r'$\delta v$ (km s$^{-1})$')

    plotcounter = plotcounter + 1
    namecounter = namecounter +1
plt.subplots_adjust(hspace=.0)
fig.text(0.02, 0.55, r'$\tilde{F}$', ha='center', va='center', rotation='vertical')
plt.show()
#
# # the fisrt subplot
# ax0 = plt.subplot(gs[0])
# #ax0.set_yscale("log")
# ax0.plot(np.array([oursvrelbins[0],oursvrelbins[-1]]),np.array([1,1]),'#FFA500', linestyle =':')
# ax0.plot(np.array([1,1]), np.array([-1,5]),'#FFA500', linestyle =':')
# ax0.step(oursvrelbins, oursmeanbinned, color='#ff0000')
# ax0.step(oursvrelbins, ourscontrolmeanbinned, color='#00ff00')
# ax0.text(24500, 1.1,'$N_{QSO}^{LOS}$ = '+str(oursqsonumber))
# ax0.text(24500, 1.05,'$N_{CARLA}$ = '+str(ourscarlanumber))
# ax0.text(24000, 0.65,'$N_{QSO}^{Control}$ = '+str(ourscontrolqsonumber))
#
# #the second subplot
# ax1 = plt.subplot(gs[2], sharex = ax0, sharey = ax0)
# ax1.plot(np.array([oursvrelbins[0],oursvrelbins[-1]]),np.array([1,1]),'#FFA500', linestyle =':')
# ax1.plot(np.array([1,1]), np.array([-1,5]),'#FFA500', linestyle =':')
# ax1.step(oursvrelbins, oursmedbinned, color='#0000ff')
# ax1.step(oursvrelbins, ourscontrolmedbinned, color='#00ff00')
# ax1.set_xlabel(r'$\delta v$ (km s$^{-1})$')
# plt.setp(ax0.get_xticklabels(), visible=False)
# # remove last tick label for the second subplot
# yticks = ax1.yaxis.get_major_ticks()
# yticks[-1].label1.set_visible(False)
#
# # the 3rd subplot
# ax2 = plt.subplot(gs[1], sharex = ax0, sharey = ax0)
# #ax0.set_yscale("log")
# ax2.plot(np.array([oursvrelbins[0],oursvrelbins[-1]]),np.array([1,1]),'#FFA500', linestyle =':')
# ax2.plot(np.array([1,1]), np.array([-1,5]),'#FFA500', linestyle =':')
# ax2.step(leesvrelbins, leesmeanbinned, color='#ff0000')
# ax2.step(leesvrelbins, leescontrolmeanbinned, color='#00ff00')
# ax2.text(24500, 1.1,'$N_{QSO}^{LOS}$ = '+str(leesqsonumber))
# ax2.text(24500, 1.05,'$N_{CARLA}$ = '+str(leescarlanumber))
# ax2.text(24000, 0.65,'$N_{QSO}^{Control}$ = '+str(leescontrolqsonumber))
#
# #the 4th subplot
# ax3 = plt.subplot(gs[3], sharex = ax0, sharey = ax0)
# ax3.plot(np.array([oursvrelbins[0],oursvrelbins[-1]]),np.array([1,1]),'#FFA500', linestyle =':')
# ax3.plot(np.array([1,1]), np.array([-1,5]),'#FFA500', linestyle =':')
# ax3.step(leesvrelbins, leesmedbinned, color='#0000ff')
# ax3.step(leesvrelbins, leescontrolmedbinned, color='#00ff00')
# ax3.set_xlabel(r'$\delta v$ (km s$^{-1})$')
# plt.setp(ax2.get_xticklabels(), visible=False)
# # remove last tick label for the second subplot
# yticks = ax3.yaxis.get_major_ticks()
# yticks[-1].label1.set_visible(False)
# ax3.set_ylim([0.6, 1.15])
# ax3.set_xlim([-40000, 60000])
#
#
# fig.text(0.05, 0.55, r'$\tilde{F}$', ha='center', va='center', rotation='vertical')
#
# # remove vertical gap between subplots
# plt.subplots_adjust(hspace=.0)
# plt.show()
