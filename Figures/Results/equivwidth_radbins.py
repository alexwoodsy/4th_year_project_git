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
#fig = plt.figure('absline fitting')
# set height ratios for sublots
#gs = plt.GridSpec(2, 3, height_ratios=[1, 1])


plotcounter = 0
namecounter = 0

runs = ['ours_200to400bin']
panelname = [r'S/N $>$ 4']
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


    ###absorbtion line fitting####
    #fit controll level
    def polynomial(x, a, b):
        return a*(x)+b

    lower, upper = -2000, 2000 #range over which the fit takes place
    plotdatarange = np.linspace(lower, upper,1000)
    contdatarange = np.arange(findval(vrelbins,lower), findval(vrelbins,upper))
    contpopt, pcov = cf(polynomial, vrelbins[contdatarange], controlmeanbinned[contdatarange], bounds =([-100,-1000],[100,1000]))
    conta, contb = contpopt
    control_fit = polynomial(vrelbins[contdatarange], *contpopt)
    control_level = 1
    #fitting needs initial data so extract data about abs line
    def guassian(x, amp, mean, std):
        return amp*np.exp(-((x-mean)**2)/(2*(std)**2))

    figabs, ax = plt.subplots(1,1,num='Absorption line fitting')
    datarange = np.arange(findval(vrelbins,lower), findval(vrelbins,upper))
    popt, pcov = cf(guassian, vrelbins[datarange], control_fit-meanbinned[datarange], bounds =([-1,-700,0],[1.5,700,1000]))
    meanamp, meanmean, meanstd = popt

    controlfitplot = polynomial(plotdatarange, *contpopt)
    guassianplot = guassian(plotdatarange, *popt)

    fittedline = controlfitplot-guassianplot


    ax.step(vrelbins, meanbinned,label='stack data')
    ax.plot(plotdatarange, fittedline , 'r-',label='fitting parmaters: amp=%5.3f, mean=%5.3f, std=%5.3f' % tuple(popt))

    ax.step(vrelbins, controlmeanbinned,label='control data')
    ax.plot(plotdatarange, controlfitplot, 'b-',label='fitting parmaters: a=%5.3f, b=%5.3f' % tuple(contpopt))

    ax.set_xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
    ax.set_ylabel(r'$<F>$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')
    ax.legend()
    ax.set_ylim([0.55, 1.05])
    ax.set_xlim([-2500, 3000])


plt.show()
