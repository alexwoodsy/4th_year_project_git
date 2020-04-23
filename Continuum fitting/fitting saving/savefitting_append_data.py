import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
import os, time
#import fitting
import sys
#sys.path.append('C:/Users/jason/GIT/4th_year_project_git/Continuum Fitting')
sys.path.append('C:/Users/alexw/Documents/GitHub/4th_year_project_git/Continuum fitting/Fitting new')
import fastfit_v1 as fitmethod

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
    ind = np.argsort(array)
    sortedarray = array[ind]
    selectedinds = ind[-pct:] #indices of the top pct in array

    start = int(len(selectedinds)/2 - pct/4)
    end = int(len(selectedinds)/2 + pct/4)
    selectedinds = selectedinds[start:end]


    #selectedinds = ind[int(-pct*0.5)]
    return selectedinds


def findpctmean(array,pct): #find median pct of data in array
#e.g if pct = 10% of array it takes 45-55 percentile
    sortedarray = np.argsort(array)
    start = int(len(sortedarray)/2 - pct/2)
    end = int(len(sortedarray)/2 + pct/2)
    selectedvals = sortedarray[start:end]
    return selectedvals


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

for i in range(3000,fitlen):
    lyalpha = metalyalpha[i]
    lyalphaind = int(metalyalphaind[i])
    prefitspec = metaspecfilename[i]
    clustername = metagcname[i]

    spec = prefitspec[0:20]+'.fits'
    specdirectory = 'Spectra/'+spec
    data = fits.getdata(specdirectory,ext=1)
    speclen = len(data)
    wlen = 10**data.field(1)
    flux = data.field(0)

    #--------------------Continuum fitting-----------------------#

    #split the spec in two about lyalpha peak for 2 fitting regions
    medintervalwlen = np.array([])
    medwinpeak = np.array([])
    maxintervalwlen = np.array([])
    maxwinpeak = np.array([])
    forestwinnum, forestpct = 20, 0.2 #forest number must be even #metadata
    otherwinnum, otherpct =  50, 0.2 #metadata

    pw = 30


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
        while step <= speclen:
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

    #pad interval with start/end value to allign correctly
    medintervalwlen[0] = wlen[0]
    medintervalwlen[-1] = wlen[-1]
    maxintervalwlen[0] = wlen[0]
    maxintervalwlen[-1] = wlen[-1]

    medintpol = interpolate.interp1d(medintervalwlen, medwinpeak, kind=1)
    maxintpol = interpolate.interp1d(maxintervalwlen, maxwinpeak, kind=1)

    medcontfit = medintpol(wlen)
    maxcontfit = maxintpol(wlen)
    #smooth fit
    if stackstatus == 'SUCCESS':
        smoothwin = 2*int(lyalphaind/forestwinnum)+1 #ensures smooth window is odd number
        if smoothwin < 10:
            smoothwin = 2*int((speclen-lyalphaind)/otherwinnum)+1
    else:
        smoothwin = 2*int((speclen-lyalphaind)/otherwinnum)+1

    medcontfit = signal.savgol_filter(medcontfit, smoothwin,3)
    maxcontfit = signal.savgol_filter(maxcontfit, smoothwin,3)


    # plt.figure()
    # plt.plot(wlen, flux)
    # plt.plot(wlen, medcontfit)
    # #plt.plot(maxintervalwlen, maxwinpeak,'*')
    # plt.plot(wlen, maxcontfit)
    # plt.plot(wlen,postlyavg)
    # plt.show()


    wlencol = fits.Column(name='Wavelength', array=wlen, format='F')
    fluxcol = fits.Column(name='Flux', array=flux, format='F')
    med_continuumcol = fits.Column(name='med_continuum', array=medcontfit, format='F')
    max_continuumcol = fits.Column(name='max_continuum', array=maxcontfit, format='F')

    specfitdata = fits.BinTableHDU.from_columns([wlencol, fluxcol, med_continuumcol, max_continuumcol])


    #get old fitsfile metadata and reappend
    hduold = fits.open('Fitted Spectra/' + prefitspec)
    oldcols = hduold[2].columns
    metadata = fits.BinTableHDU.from_columns(oldcols)
    hduold.close()

    primary = fits.PrimaryHDU()
    hdul = fits.HDUList([primary, specfitdata, metadata])

    outname = 'Fitted Spectra/' + prefitspec

    hdul.writeto(outname, overwrite = True)
    print('done '+prefitspec + ' number ' + str(i))








#
