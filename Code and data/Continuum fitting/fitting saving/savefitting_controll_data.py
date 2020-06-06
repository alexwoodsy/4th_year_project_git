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

specnames = next(os.walk('E:/Sky Spectra Zip/Sky Spectra'))[2]
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
#get precalculate metadata needed for fitting and getting spec:
metadata = fits.getdata('Anti-Match/Anti-Match-Metadata.fits',ext=1)
metaplate = metadata.field(5)
metamjd = metadata.field(6)
metafiberid = metadata.field(7)
metaredshift = metadata.field(9)

#read in spec and do the fitting
for i in range(20000,len(metaplate)):
    plate = str(metaplate[i]).zfill(4)
    mjd = str(metamjd[i]).zfill(5)
    fiberid = str(metafiberid[i]).zfill(4)
    spec = 'spec-'+plate+'-'+mjd+'-'+fiberid+'.fits'
    specdirectory = 'E:/Sky Spectra Zip/Sky Spectra/'+spec

#check spec has been downloaded
    if os.path.exists(specdirectory) == False:
        print(spec+' was not downloaded!')
    else:

        data = fits.getdata(specdirectory,ext=1)
        speclen = len(data)
        wlen = 10**data.field(1)
        flux = data.field(0)

        #get meta inof needed for cont fit for this stack
        redshift = metaredshift[i]
        lyalpha = (1+redshift)*1215.67
        lyalphaind = findval(wlen, lyalpha)




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

        #saving
        wlencol = fits.Column(name='Wavelength', array=wlen, format='F')
        fluxcol = fits.Column(name='Flux', array=flux, format='F')
        # med_continuumcol = fits.Column(name='med_continuum_2', array=medcontfit, format='F')
        # max_continuumcol = fits.Column(name='max_continuum_2', array=maxcontfit, format='F')
        max_medcontinuumcol = fits.Column(name='med_max_continuum', array=max_medcontfit, format='F')

        specfitdata = fits.BinTableHDU.from_columns([wlencol, fluxcol, max_medcontinuumcol])


        primary = fits.PrimaryHDU()
        hdul = fits.HDUList([primary, specfitdata])
        name = 'AM-spec-'+plate+'-'+mjd+'-'+fiberid+'-Preffited.fits'
        outname = 'E:/AM-Prefitted Spectra/' + name

        hdul.writeto(outname, overwrite = True)
        print('done '+ name + ' number ' + str(i))








#
