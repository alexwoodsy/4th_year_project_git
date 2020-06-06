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

for i in range(3000,fitlen):
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



    #
    #
    #
    # plt.figure('med cont')
    # plt.plot(wlen, flux)
    # plt.plot(wlen, medcontfit)
    # plt.plot(medintervalwlen, medwinpeak,'*')
    #
    # plt.figure('max cont')
    # plt.plot(wlen, flux)
    # plt.plot(wlen, maxcontfit)
    # plt.plot(maxintervalwlen, maxwinpeak,'*')
    # plt.plot(wlenline, 10+fluxline,'.')
    #
    # plt.figure('max_med cont')
    # plt.plot(wlen, flux)
    # plt.plot(wlen, max_medcontfit)
    # plt.plot(max_medintervalwlen, max_medwinpeak,'*')
    # plt.plot(wlenline, 10+fluxline,'.')
    #
    #
    # plt.show()


    #old data in
    hduold = fits.open('Fitted Spectra/' + prefitspec)
    specfitoldcols = hduold[1].columns

    #bin our data and append it to the old fitting we did
    med_continuumcol = fits.Column(name='med_continuum_2', array=medcontfit, format='F')
    max_continuumcol = fits.Column(name='max_continuum_2', array=maxcontfit, format='F')
    max_medcontinuumcol = fits.Column(name='med_max_continuum', array=max_medcontfit, format='F')

    specfitdata = fits.BinTableHDU.from_columns(specfitoldcols+med_continuumcol+max_continuumcol+max_medcontinuumcol)

    hduold.close()


    #read in lees continuum fit if possible:
    #get lees fitting data:
    plate = spec[5:9]
    mjd = spec[10:15]
    fiberid = spec[16:20]

    leepath = 'E:/spectralyalpha/BOSSLyaDR9_spectra/BOSSLyaDR9_spectra/'+plate+'/'+'speclya-'+plate+'-'+mjd+'-'+fiberid+'.fits'


    if os.path.exists(leepath) == False: #check lee fit is in sample, if not output just ours +flags
        #get old fitsfile metadata and reappend
        hduold = fits.open('Fitted Spectra/' + prefitspec)
        metaoldcols = hduold[2].columns
        #SDSS NAME
        leecheckcol = fits.Column(name='LEE_CHECK', array = np.array([0]), format='K')
        leecontflagcol = fits.Column(name='LEE_CONT_FLAG', array = np.array([0]), format='K')

        metadata = fits.BinTableHDU.from_columns(metaoldcols+leecheckcol+leecontflagcol)
        hduold.close()

        primary = fits.PrimaryHDU()
        hdul = fits.HDUList([primary, specfitdata, metadata])

        outname = 'Fitted Spectra/' + prefitspec

        hdul.writeto(outname, overwrite = True)
        print('done '+prefitspec + ' number {lee}' + str(i))
    else:
        #get lee continuum fit
        hdul = fits.open(leepath)
        leedata = hdul[1].data
        contflag = hdul[1].header['CONTFLG']

        #leedata = fits.getdata(path,ext=1)#import fits image
        leeflux = leedata.field(0)
        leewlen = 10**leedata.field(1)
        leemodel = leedata.field(7)
        noisecorr = leedata.field(10)
        leecont = leedata.field(11)

        hdul.close()

        #BIN LEE DATA
        leewlencol = fits.Column(name='Wavelength', array=leewlen, format='F')
        leefluxcol = fits.Column(name='Flux', array=leeflux, format='F')
        lee_contcol = fits.Column(name='Continuum', array=leecont, format='F')
        lee_noisecorrcol = fits.Column(name='Noise_Correction', array=noisecorr, format='F')

        leefitdata = fits.BinTableHDU.from_columns([leewlencol, leefluxcol, lee_contcol, lee_noisecorrcol])

        #get old fitsfile metadata and reappend
        hduold = fits.open('Fitted Spectra/' + prefitspec)
        metaoldcols = hduold[2].columns
        #SDSS NAME
        leecheckcol = fits.Column(name='LEE_CHECK', array = np.array([1]).astype(int), format='K')
        leecontflagcol = fits.Column(name='LEE_CONT_FLAG', array = np.array([contflag]).astype(int), format='K')
        metadata = fits.BinTableHDU.from_columns(metaoldcols+leecheckcol+leecontflagcol)
        hduold.close()

        primary = fits.PrimaryHDU()
        hdul = fits.HDUList([primary, specfitdata, leefitdata, metadata])

        outname = 'Fitted Spectra/' + prefitspec

        hdul.writeto(outname, overwrite = True)
        print('done '+prefitspec + ' number ' + str(i))








#
