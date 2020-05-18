import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
import os, time
#import fitting
import sys
#sys.path.append('C:/Users/jason/GIT/4th_year_project_git/Continuum Fitting')
plt.style.use('mystyle')

specnames = next(os.walk('Spectra/'))[2]
spectot = len(specnames)



#definitions for use in fitting
def findval(array,val): # find index of nearst value in array
    array = np.asarray(array)
    ind = np.abs(array - val).argmin()
    return ind

def findmax(array): #find max index
    array = np.asarray(array)
    ind = (np.abs(array - np.max(array))).argmin()
    return ind

def findmed(array): #find median index
    array = np.asarray(array)
    ind = (np.abs(array - np.median(array))).argmin()
    return ind

def findpctmax(array,pct):
    sortedarray = np.argsort(array)
    if len(sortedarray) < pct:
        selectedinds = sortedarray[int(len(sortedarray)*0.9)]
    else:
        selectedinds = sortedarray[-pct:] #indices of the top pct in array
        if len(selectedinds) != 0:
            selectedinds = selectedinds[int(-pct*0.9)] #take single value

    # start = int(len(selectedinds)/2 - pct/4)
    # end = int(len(selectedinds)/2 + pct/4)
    # selectedinds = selectedinds[start:end]

    return selectedinds


def findpctmean(array,pct): #find median pct of data in array
#e.g if pct = 10% of array it takes 45-55 percentile
    sortedarray = np.argsort(array)
    start = int(len(sortedarray)/2 - pct/2)
    end = int(len(sortedarray)/2 + pct/2)
    selectedinds = sortedarray[start:end]
    if len(selectedinds) != 0:
        selectedinds = sortedarray[int(len(sortedarray)/2)]

    return selectedinds


#---------------------------data extraction-----------------------------#
#get precalculate metadata needed for fitting:
#import metadata tableand data
metatabledata = fits.getdata('metatable.fits', ext=1)
fitlen = len(metatabledata)

metaspecfilename = []
metaredshift = np.zeros(fitlen)
metastonall = np.zeros(fitlen)
metastonforest = np.zeros(fitlen)
metalyalpha = np.zeros(fitlen)
metalyalphaind = np.zeros(fitlen)
metagcname = []
metagcredshift = np.zeros(fitlen)
metagc_qso_sep = np.zeros(fitlen)
metagclyalpha = np.zeros(fitlen)
metagclyalphaind = np.zeros(fitlen)
metastackmsg = []


for i in range(0,fitlen):
    #extract metadata
    metaspecfilename.append(metatabledata[i][0])
    metaredshift[i]= metatabledata[i][1]
    metastonall[i] = metatabledata[i][2]
    metastonforest[i] = metatabledata[i][3]
    metalyalpha[i] = metatabledata[i][4]
    metalyalphaind[i] = metatabledata[i][5]
    metagcname.append(metatabledata[i][6])
    metagcredshift[i] = metatabledata[i][7]
    metagc_qso_sep[i] = metatabledata[i][8]
    metagclyalpha[i] = metatabledata[i][9]
    metagclyalphaind[i] = metatabledata[i][10]
    metastackmsg.append(metatabledata[i][11])



fig = plt.figure('fitting')
# set height ratios for sublots
gs = plt.GridSpec(3, 2, height_ratios=[1, 1, 1])
plotcounter = 0

for i in [0,10,1000]: #+10,1000
    lyalpha = metalyalpha[i]
    redshift = metaredshift[i]
    lyalphaind = int(metalyalphaind[i])
    prefitspec = metaspecfilename[i]
    clustername = metagcname[i]

    spec = prefitspec[0:20]+'.fits'
    specdirectory = 'Spectra/'+spec
    data = fits.getdata(specdirectory,ext=1)
    speclen = len(data)
    wlen = 10**data.field(1)
    flux = data.field(0)


    #get spec name from old data:
    sdssdata = fits.getdata('PositionsTable.fits',ext=1)
    i = 0
    namecheck = False
    for i in range(len(sdssdata)):
        file = 'spec-'+str(sdssdata[i][4]).zfill(4) +'-'+ str(sdssdata[i][5]).zfill(5)+'-' + str(sdssdata[i][6]).zfill(4)
        if str(file) == spec[0:20]:
            sdssspecname = sdssdata[i][0]



    #--------------------Continuum fitting-----------------------#

    #split the spec in two about lyalpha peak for 2 fitting regions
    medintervalwlen = np.array([])
    medwinpeak = np.array([])
    maxintervalwlen = np.array([])
    maxwinpeak = np.array([])
    forestwinnum, forestpct = 20, 0.4 #forest number must be even #metadata
    otherwinnum, otherpct =  50, 0.4 #metadata

    pw = 100

    #fitv9 method
    if lyalpha - wlen[0] <= pw: #if no forest just fit the other part of the spectrum
        stackstatus = 'FORESTERROR' #metadata
        step = 0
        while step <= speclen:
            window = int((speclen-lyalphaind)/otherwinnum)
            percentage = otherpct
            pct = int(window*percentage)
            windata = flux[step:(step+window)]

            #med
            medwinpeakind = step + findpctmean(windata,pct)
            medwinpeak = np.append(medwinpeak,flux[medwinpeakind])
            medintervalwlen = np.append(medintervalwlen,wlen[medwinpeakind])
            #max
            maxwinpeakind = step + findpctmax(windata,pct)
            maxwinpeak = np.append(maxwinpeak,flux[maxwinpeakind])
            maxintervalwlen = np.append(maxintervalwlen,wlen[maxwinpeakind])

            step = step + window
    else: #for good spec fit accordingly
        stackstatus = 'SUCCESS' #metadata
        step = 0
        while step < speclen:
            if step <= lyalphaind:
                winnum = forestwinnum
                window = int(lyalphaind/winnum)
                # if window < lyalphaind:
                #     window = int((speclen-lyalphaind)/otherwinnum)
                percentage = forestpct
                pct = int(window*percentage)
                windata = flux[step:(step+window)]

                #med
                medwinpeakind = step + findpctmean(windata,pct)
                medwinpeak = np.append(medwinpeak,flux[medwinpeakind])
                medintervalwlen = np.append(medintervalwlen,wlen[medwinpeakind])
                #max
                maxwinpeakind = step + findpctmax(windata,pct)
                maxwinpeak = np.append(maxwinpeak,flux[maxwinpeakind])
                maxintervalwlen = np.append(maxintervalwlen,wlen[maxwinpeakind])

                step = step + window
            else:
                winnum = otherwinnum
                window = int((speclen-lyalphaind)/winnum)
                percentage = otherpct
                pct = int(window*percentage)
                windata = flux[step:(step+window)]

                #med
                medwinpeakind = step + findpctmean(windata,pct)
                medwinpeak = np.append(medwinpeak,flux[medwinpeakind])
                medintervalwlen = np.append(medintervalwlen,wlen[medwinpeakind])
                #max
                maxwinpeakind = step + findpctmax(windata,pct)
                maxwinpeak = np.append(maxwinpeak,flux[maxwinpeakind])
                maxintervalwlen = np.append(maxintervalwlen,wlen[maxwinpeakind])

                step = step + window


    #add in points for splines at emmission lines
    refline = open('refline.txt', 'r')
    linename = []
    wlenline = np.array([])
    for line in refline:
        line = line.strip()
        columns = line.split()
        wlenline = np.append(wlenline, float(columns[0]))
        linename.append(columns[1])
    refline.close()
    linename = linename[0:12] #restrict higher WLEN lines and 1st line out of range
    wlenline = wlenline[0:12]*(redshift+1)
    fluxline = np.zeros(len(wlenline))
    for j in range(0,len(wlenline)):
        fluxline[j] = flux[findval(wlenline[j],wlen)]

    #remove point too close
    for peak in wlenline:
        medpeakdelinds = []
        maxpeakdelinds = []
        for ind in range(0,len(maxintervalwlen)):
            if abs(peak - maxintervalwlen[ind]) <= 10:
                maxpeakdelinds.append(ind)
        maxintervalwlen = np.delete(maxintervalwlen, maxpeakdelinds)
        maxwinpeak = np.delete(maxwinpeak, maxpeakdelinds)

        for ind in range(0,len(medintervalwlen)):
            if abs(peak - medintervalwlen[ind]) <= 10:
                medpeakdelinds.append(ind)
        medintervalwlen = np.delete(medintervalwlen, medpeakdelinds)
        medwinpeak = np.delete(medwinpeak, medpeakdelinds)

    #add in emmision peaks
    maxintervalwlen = np.append(maxintervalwlen,wlenline)
    sortind = np.argsort(maxintervalwlen)
    maxintervalwlen = maxintervalwlen[sortind]
    maxwinpeak = np.append(maxwinpeak, fluxline)
    maxwinpeak = maxwinpeak[sortind]

    medintervalwlen = np.append(medintervalwlen,wlenline)
    sortind = np.argsort(medintervalwlen)
    medintervalwlen = medintervalwlen[sortind]
    medwinpeak = np.append(medwinpeak, fluxline)
    medwinpeak = medwinpeak[sortind]

    #pad interval with start/end value to allign correctly
    medintervalwlen[0] = wlen[0]
    medintervalwlen[-1] = wlen[-1]
    maxintervalwlen[0] = wlen[0]
    maxintervalwlen[-1] = wlen[-1]

    #combined where max is used in the forest up to lyalpha peak and med after;
    maxlyaintind = findval(maxintervalwlen, wlenline[1])
    medlyaintind = findval(maxintervalwlen, wlenline[1])
    max_medintervalwlen = np.array([])
    max_medintervalwlen = np.append(maxintervalwlen[0:maxlyaintind],medintervalwlen[medlyaintind:])
    max_medwinpeak = np.append(maxwinpeak[0:maxlyaintind],medwinpeak[medlyaintind:])

    # plt.plot(wlen,flux)
    # plt.plot(medintervalwlen, medwinpeak)
    # plt.show()


    medintpol = interpolate.interp1d(medintervalwlen, medwinpeak, kind = 1)

    maxintpol = interpolate.interp1d(maxintervalwlen, maxwinpeak, kind = 1)
    max_medintpol = interpolate.interp1d(max_medintervalwlen, max_medwinpeak, kind = 1)


    medcontfit = medintpol(wlen)
    maxcontfit = maxintpol(wlen)
    max_medcontfit = max_medintpol(wlen)

    mednorm = flux/medcontfit
    maxnorm = flux/maxcontfit
    max_mednorm = flux/max_medcontfit

    #smooth fit
    if stackstatus == 'SUCCESS':
        smoothwin = 2*int(lyalphaind/forestwinnum)+1 #ensures smooth window is odd number
        if smoothwin < 10:
            smoothwin = 2*int((speclen-lyalphaind)/otherwinnum)+1
    else:
        smoothwin = 2*int((speclen-lyalphaind)/otherwinnum)+1

    medcontfit = signal.savgol_filter(medcontfit, smoothwin,3)
    maxcontfit = signal.savgol_filter(maxcontfit, smoothwin,3)
    max_medcontfit = signal.savgol_filter(max_medcontfit, smoothwin,3)

    #topleft
    if plotcounter == 0:
        ax0 = plt.subplot(gs[plotcounter])
        ax0.plot(wlen, flux)
        ax0.plot(wlen, medcontfit)
        ax0.plot(medintervalwlen, medwinpeak,'^',color='orange')
        ax0.text(5250,135,sdssspecname)
        ax0.text(5250,115,r'z = '+str(np.round(redshift,2)))
        #ax0.plot(wlenline, fluxline,'.')

        axr = plt.subplot(gs[plotcounter+1],sharex = ax0, sharey = ax0)
        axr.plot(wlen, flux)
        axr.plot(wlen, maxcontfit)
        axr.plot(maxintervalwlen, maxwinpeak,'^')
        axr.text(5300,135,sdssspecname)
        axr.text(5300,115,r'z = '+str(np.round(redshift,2)))
        #axr.plot(wlenline, fluxline,'.')
    else:
        axl = plt.subplot(gs[plotcounter],sharex = ax0)
        axl.plot(wlen, flux)
        axl.plot(wlen, medcontfit)
        axl.plot(medintervalwlen, medwinpeak,'^',color='orange')
        #axl.plot(wlenline, fluxline,'.')

        axr = plt.subplot(gs[plotcounter+1],sharex = ax0, sharey = axl)
        axr.plot(wlen, flux)
        axr.plot(wlen, maxcontfit)
        axr.plot(maxintervalwlen, maxwinpeak,'^')
        #axr.plot(wlenline, fluxline,'.')

        if plotcounter == 2:
            axl.text(3850,12,sdssspecname)
            axl.text(3850,10,r'z = '+str(np.round(redshift,2)))
            axr.text(3850,12,sdssspecname)
            axr.text(3850,10,r'z = '+str(np.round(redshift,2)))

            axl.set_ylabel(r'$F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')

        else:
            axl.text(5250,25,sdssspecname)
            axl.text(5250,21,r'z = '+str(np.round(redshift,2)))
            axr.text(5250,25,sdssspecname)
            axr.text(5250,21,r'z = '+str(np.round(redshift,2)))

            axl.set_xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
            axr.set_xlabel(r'$\lambda$ ($\mathrm{\AA}$)')



    plotcounter = plotcounter + 2

# remove vertical gap between subplots
plt.subplots_adjust(hspace=.0)

#labels
# Set common labels
#fig.text(0.555, 0.03, r'$\lambda$ ($\mathrm{\AA}$)', ha='center', va='center')
#right
#fig.text(0.58, 0.5, r'$\tilde{F}$', ha='center', va='center', rotation='vertical')
#left
#fig.text(0.05, 0.5, r'$F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$', ha='center', va='center', rotation='vertical')
plt.show()

    #For plotting all in 1 go
    # fig = plt.figure('fitting')
    # # set height ratios for sublots
    # gs = plt.GridSpec(3, 2, height_ratios=[1, 1, 1])
    # # the fisrt subplot
    # ax00 = plt.subplot(gs[0])
    # #ax0.set_yscale("log")
    # ax00.plot(wlen, flux)
    # ax00.plot(wlen, maxcontfit)
    # ax00.plot(maxintervalwlen, maxwinpeak,'*')
    # ax00.plot(wlenline, fluxline,'.')
    # #name0 = ax0.text(8750,140,sdssspecname)
    # #z0 = ax0.text(9500,120,r'z = '+str(z))
    # ax01 = plt.subplot(gs[2], sharex = ax00)
    # ax01.plot(wlen, flux)
    # ax01.plot(wlen, medcontfit)
    # ax01.plot(medintervalwlen, medwinpeak,'*')
    # ax01.plot(wlenline, fluxline,'.')
    #
    # ax02 = plt.subplot(gs[4], sharex = ax00)
    # ax02.plot(wlen, flux)
    # ax02.plot(wlen, max_medcontfit)
    # ax02.plot(max_medintervalwlen, max_medwinpeak,'*')
    # ax02.plot(wlenline, fluxline,'.')
    #
    # # the fisrt subplot
    # ax10 = plt.subplot(gs[1], sharex = ax00)
    # #ax0.set_yscale("log")
    # ax10.plot(wlen, maxnorm)
    #
    #
    # #name0 = ax0.text(8750,140,sdssspecname)
    # #z0 = ax0.text(9500,120,r'z = '+str(z))
    # ax11 = plt.subplot(gs[3], sharex = ax10)
    # ax11.plot(wlen, mednorm)
    #
    #
    # ax12 = plt.subplot(gs[5], sharex = ax10)
    # ax12.plot(wlen, max_mednorm)
    #
    #
    # # remove vertical gap between subplots
    # plt.subplots_adjust(hspace=.0)
    # #labels
    # # Set common labels
    # fig.text(0.5, 0.03, r'$\lambda$ ($\mathrm{\AA}$)', ha='center', va='center')
    # #right
    # fig.text(0.58, 0.5, r'$\tilde{F}$', ha='center', va='center', rotation='vertical')
    # #left
    # fig.text(0.05, 0.5, r'$F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$', ha='center', va='center', rotation='vertical')
    # plt.show()
