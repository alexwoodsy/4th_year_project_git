import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal, integrate
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
gs = plt.GridSpec(1, 3, height_ratios=[1])


plotcounter = 0
namecounter = 0

runs = ['ours_200to400bin_sn_2','ours_200to400bin_odl2_sn2','ours_200to400bin_odg2_sn2']
panelname = [r'all $\sigma$',r'$\sigma <$ 2',r'$\sigma \geq$ 2']
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
    #print('qso number  = '+str(qsonumber))
    #print('carla number  = '+str(carlanumber))
    #print('control number  = '+str(controlqsonumber))

    step =0
    for i in range(0,len(vrelbins)):
        binind = findval(vrel, vrelbins[i])
        meanbinned[i] = np.mean(mean[step:binind])
        medbinned[i] = np.median(med[step:binind])
        controlmeanbinned[i] = np.mean(controlmean[step:binind])
        controlmedbinned[i] = np.median(controlmed[step:binind])
        step = binind


    ###absorbtion line fitting####
    #fit controll level
    def polynomial(x, a, b):
        return a*(x)+b



    lower, upper = -4000, 4000 #range over which the fit takes place

    plotdatarange = np.linspace(lower, upper,1000)
    contdatarange = np.arange(findval(vrelbins,lower), findval(vrelbins,upper))
    contpopt, pcov = cf(polynomial, vrelbins[contdatarange], controlmeanbinned[contdatarange], bounds =([-100,-1000],[100,1000]))
    conta, contb = contpopt
    control_fit = polynomial(vrelbins[contdatarange], *contpopt)
    control_level = 1
    #fitting needs initial data so extract data about abs line
    def guassian(x, amp, mean, std):
        return amp*np.exp(-((x-mean)**2)/(2*(std)**2))


    datarange = np.arange(findval(vrelbins,lower), findval(vrelbins,upper))
    popt, pcov = cf(guassian, vrelbins[datarange], control_fit-meanbinned[datarange], bounds =([0,-1500,0],[1.5,1500,1000]))

    meanamp, meanmean, meanstd = popt
    controlfitplot = polynomial(plotdatarange, *contpopt)
    guassianplot = guassian(plotdatarange, *popt)

    fittedline = controlfitplot-guassianplot

    #convert back to wavelength:
    wlenrange = (plotdatarange/c)*lam_em + lam_em
    area = integrate.simps(guassianplot, wlenrange)
    width = area/meanamp
    col_density = (1.84e14)*width
    print('fitting parmaters: amp=%5.3f, mean=%5.3f, std=%5.3f' % tuple(popt))
    print('equiv width = '+str(width))
    print('column density = '+"{:.2e}".format(col_density))


    # the fisrt subplot seperately
    if plotcounter == 0:
        axr0 = plt.subplot(gs[plotcounter], label=runsavename)

        axr0.step((vrelbins+(binsize/2)), meanbinned, color='#0000ff', label='Composite flux')
        axr0.plot(plotdatarange, fittedline , '#ff0000',label='Gaussian fit')

        axr0.step((vrelbins+(binsize/2)), controlmeanbinned, color='#00ff00',label='Control')
        axr0.plot(plotdatarange, controlfitplot, '#FF1493', linestyle =':',linewidth=1.5,label='Control level fit')

        axr0.text(500, 0.6,panelname[namecounter])
        axr0.set_ylim([0.5, 1.05])
        axr0.set_xlim([-700, 700])
        axr0.set_xlabel(r'$\delta v$ (km s$^{-1})$')
    else:

        axr = plt.subplot(gs[plotcounter],sharex = axr0, sharey = axr0, label=runsavename)

        axr.step((vrelbins+(binsize/2)), meanbinned, color='#0000ff',label='Composite flux')
        axr.plot(plotdatarange, fittedline , '#ff0000',label='Gaussian fit')

        axr.step((vrelbins+(binsize/2)), controlmeanbinned, color='#00ff00',label='Control')
        axr.plot(plotdatarange, controlfitplot, '#FF1493', linestyle =':',linewidth=1.5,label='Control level fit')
        axr.text(500, 0.6,panelname[namecounter])

        axr.set_ylim([0.5, 1.05])
        axr.set_xlim([-700, 700])
        if runsavename != runs[-1]: #
            yticks = axr.yaxis.get_major_ticks()
            xticks = axr.xaxis.get_major_ticks()
            yticks[-1].label1.set_visible(False)
            xticks[-1].label1.set_visible(False)
            axr.set_xlabel(r'$\delta v$ (km s$^{-1})$')
        else:
            axr.set_xlabel(r'$\delta v$ (km s$^{-1})$')


    #
    if namecounter == 0:
        axr0.set_ylabel(r'$\tilde{F}$')
        leg = plt.legend()
        leg.get_frame().set_linewidth(0.0)

    plotcounter = plotcounter + 1
    namecounter = namecounter +1
plt.subplots_adjust(hspace=.0)
#fig.text(0.05, 0.55, r'$\tilde{F}$', ha='center', va='center', rotation='vertical')
plt.show()

#
#     ###absorbtion line fitting####
#     #fit controll level
#     def polynomial(x, a, b):
#         return a*(x)+b
#
#     lower, upper = -4000, 4000 #range over which the fit takes place
#     plotdatarange = np.linspace(lower, upper,1000)
#     contdatarange = np.arange(findval(vrelbins,lower), findval(vrelbins,upper))
#     contpopt, pcov = cf(polynomial, vrelbins[contdatarange], controlmeanbinned[contdatarange], bounds =([-100,-1000],[100,1000]))
#     conta, contb = contpopt
#     control_fit = polynomial(vrelbins[contdatarange], *contpopt)
#     control_level = 1
#     #fitting needs initial data so extract data about abs line
#     def guassian(x, amp, mean, std):
#         return amp*np.exp(-((x-mean)**2)/(2*(std)**2))
#
#     figabs, ax = plt.subplots(1,1,num='Absorption line fitting')
#     datarange = np.arange(findval(vrelbins,lower), findval(vrelbins,upper))
#     popt, pcov = cf(guassian, vrelbins[datarange], control_fit-meanbinned[datarange], bounds =([-1,-1500,0],[1.5,1500,1000]))
#     meanamp, meanmean, meanstd = popt
#     controlfitplot = polynomial(plotdatarange, *contpopt)
#     guassianplot = guassian(plotdatarange, *popt)
#
#     fittedline = controlfitplot-guassianplot
#     wlenbins = ((vrelbins+(binsize/2))/c)*lam_em + lam_em
#     wlenplotrange = (plotdatarange/c)*lam_em + lam_em
#
#     ax.step((vrelbins+(binsize/2)), meanbinned, color='#0000ff', label=panelname[namecounter])
#     ax.plot(plotdatarange, fittedline , '#ff0000',label='Gaussian fit')
#
#     ax.step((vrelbins+(binsize/2)), controlmeanbinned, color='#00ff00',label='Control')
#     ax.plot(plotdatarange, controlfitplot, '#FF1493', linestyle =':',linewidth=1.5,label='Control level fit')
#
#     ax.set_xlabel(r'$\delta v$ (km s$^{-1})$')
#     ax.set_ylabel(r'$\tilde{F}$')
#     leg = plt.legend()
#     leg.get_frame().set_linewidth(0.0)
#     ax.set_ylim([0.55, 1.05])
#     ax.set_xlim([-2500, 3000])
#     #convert back to wavelength:
#     wlenrange = (plotdatarange/c)*lam_em + lam_em
#     area = integrate.simps(guassianplot, wlenrange)
#     width = area/meanamp
#     col_density = (1.84e14)*width
#     print('fitting parmaters: amp=%5.3f, mean=%5.3f, std=%5.3f' % tuple(popt))
#     print('equiv width = '+str(width))
#     print('dolumn density = '+"{:.2e}".format(col_density))
#
#     #plt.figure('test')
#     #plt.plot(wlenrange,guassianplot)
#     namecounter = namecounter + 1
# plt.show()
