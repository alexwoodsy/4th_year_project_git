import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
import os
#sys.path.append('C:/Users/alexw/Documents/GitHub/4th_year_project_git/Continuum fitting/Fitting new')
import fastfit_v1 as fitmethod

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

def findpctmax(array,pct): #find 50th percentile of top pct % in array
    ind = np.argsort(array)
    sortedarray = array[ind]
    start = int(len(sortedarray)/2 - pct/2)
    end = int(len(sortedarray)/2 + pct/2)
    selectedvals = sortedarray[start:end]
    #selectedinds = ind[-pct:] #indices of the top pct in array

    selectedinds = ind[int(-pct*0.5)]
    return selectedinds

def findpctmean(array,pct): #find median pct of data in array
#e.g if pct = 10% of array it takes 45-55 percentile
    sortedarray = np.argsort(array)
    start = int(len(sortedarray)/2 - pct/2)
    end = int(len(sortedarray)/2 + pct/2)
    selectedvals = sortedarray[start:end]
    return selectedvals



def fitmed(wlen, flux, lyalpha, lyalphaind):
#--------------------Continuum fitting-----------------------#

    #split the spec in two about lyalpha peak for 2 fitting regions
    speclen = len(wlen)
    pw = 30
    intervalwlen = np.array([])
    winpeak = np.array([])
    forestwinnum, forestpct = 8, 0.2 #forest number must be even #metadata
    otherwinnum, otherpct =  50, 0.2 #metadata

    #fitv9 method
    if lyalpha - wlen[0] <= 0: #if no forest just fit the other part of the spectrum
        stackstatus = 'FORESTERROR' #metadata
        step = 0
        while step <= speclen:
            window = int((speclen-lyalphaind)/otherwinnum)
            percentage = otherpct
            pct = int(window*percentage)
            windata = flux[step:(step+window)]
            winpeakind = step + findpctmean(windata,pct)
            winpeak = np.append(winpeak,flux[winpeakind])
            intervalwlen = np.append(intervalwlen,wlen[winpeakind])
            step = step + window
    else: #for good spec fit accordingly
        stackstatus = 'SUCCESS' #metadata
        step = 0
        while step <= speclen:
            if step <= lyalphaind :
                winnum = forestwinnum
                window = int(lyalphaind/winnum)
                # if window < lyalphaind:
                #     window = int((speclen-lyalphaind)/otherwinnum)
                percentage = forestpct
                pct = int(window*percentage)
                windata = flux[step:(step+window)]
                winpeakind = step + findpctmean(windata,pct)#change
                winpeak = np.append(winpeak,flux[winpeakind])
                intervalwlen = np.append(intervalwlen,wlen[winpeakind])
                step = step + window
            else:
                winnum = otherwinnum
                window = int((speclen-lyalphaind)/winnum)
                percentage = otherpct
                pct = int(window*percentage)
                windata = flux[step:(step+window)]
                winpeakind = step + findpctmean(windata,pct)
                winpeak = np.append(winpeak,flux[winpeakind])
                intervalwlen = np.append(intervalwlen,wlen[winpeakind])
                step = step + window

        #pad interval with start/end value to allign correctly
        intervalwlen[0] = wlen[0]
        intervalwlen[-1] = wlen[-1]
        intpol = interpolate.interp1d(intervalwlen, winpeak, kind=1)
        contfit = intpol(wlen)
        #smooth fit
        if stackstatus == 'SUCCESS':
            smoothwin = 2*int(lyalphaind/forestwinnum)+1 #ensures smooth window is odd number
            if smoothwin < 10:
                smoothwin = 2*int((speclen-lyalphaind)/otherwinnum)+1
        else:
            smoothwin = 2*int((speclen-lyalphaind)/otherwinnum)+1


        contfit = signal.savgol_filter(contfit, smoothwin,3)
        return contfit
