#top x% method for fitting continuum in the forest
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate
import os

#plt.style.use('mystyle') #path C:\Users\alexw\AppData\Local\Programs\Python\Python37\Lib\site-packages\matplotlib\mpl-data\stylelib

#calculates points slected in interval
def findmax(array):
    array = np.asarray(array)
    ind = (np.abs(array - np.max(array))).argmin()
    return ind

def findmed(array):
    array = np.asarray(array)
    ind = (np.abs(array - np.median(array))).argmin()
    return ind

def findpctmax(array,pct):
    sortedarray = np.argsort(array)
    selectedvals = sortedarray[-pct:]
    return selectedvals

def findpctmean(array,pct):
    sortedarray = np.argsort(array)
    start = int(len(sortedarray)/2 - pct/2)
    end = int(len(sortedarray)/2 + pct/2)
    selectedvals = sortedarray[start:end]
    return selectedvals


#imports the spectra from the spectra folder
specnames = next(os.walk('Spectra'))[2]
spectot = len(specnames)
#add indexing for spectra in file to allow loop over all

def contfitv7(specsample,zlim,stonlim,showplot,showerror):

    for spec in specsample:
        specdirectory = 'Spectra/'+spec
        #print(specdirectory)
        data = fits.getdata(specdirectory,ext=1)#import fits image
        speclen = len(data)
        flux = np.zeros(speclen)
        wlen = np.zeros(speclen)
        model = np.zeros(speclen)
        ivar = np.zeros(speclen)

        for i in range(0,speclen):
         flux[i] = data[i][0]
         ivar[i] = data[i][2]
         if ivar[i] == 0:
            ivar[i] = 0.00001
         wlen[i] = 10**(data[i][1])
         model[i] = data[i][7]


    #meta data extraction to get z:
        fitdata = fits.getdata(specdirectory,ext=2)#import fits image
        metasize = len(fitdata[0])
        #print(metasize)
        if metasize == 126:
             redshift = fitdata[0][63]
        else:
             redshift = fitdata[0][38]

    #s/n checking
        std = (1/ivar)**0.5
        ston = np.median(flux/std)


        #lyalpha considerations
        lyalphacalc = 1215.67*(1+redshift)#calc lya using redshift
        #print('lyalpha = ',lyalphacalc)
        lyalphaind = (np.abs(wlen - lyalphacalc)).argmin()#finds index of nearest point in data

    #fitting:
        #split the spec in two about lyalpha peak
        pw = 150
        intervalwlen = np.array([0])
        winpeak = np.array([0])

        #s/n check and
        #loop increments

        if redshift < zlim and showerror:
            normspec = np.zeros(speclen)
            if showerror == True:
                print('z warning!: '+spec+' z = '+ str(redshift) +' too low - normspec = 0 array')
        elif ston < stonlim:
            normspec = np.zeros(speclen)
            if showerror == True:
                print('S/N warning!: '+spec+' S/N = '+ str(ston) +' too low - normspec = 0 array')
        else:
            step = 0
            while step <= speclen:
                if step <= lyalphaind:
                    winnum = 50
                    window = int(speclen/winnum)
                    percentage = 0.2
                    pct = int(window*percentage)

                    windata = flux[step:(step+window)]
                    winpeakind = step + findpctmax(windata,pct)
                    winpeak = np.append(winpeak,flux[winpeakind])
                    intervalwlen = np.append(intervalwlen,wlen[winpeakind])
                    step = step + window
                    if lyalphaind-step < int(pw/2):
                        winpeakind = np.arange(step,(step+window))
                        winpeak = np.append(winpeak,flux[winpeakind])
                        intervalwlen = np.append(intervalwlen,wlen[winpeakind])
                        step = step + window
                else:
                    winnum = 90
                    window = int(speclen/winnum)
                    percentage = 0.2
                    pct = int(window*percentage)

                    windata = flux[step:(step+window)]
                    winpeakind = step + findpctmean(windata,pct)
                    winpeak = np.append(winpeak,flux[winpeakind])
                    intervalwlen = np.append(intervalwlen,wlen[winpeakind])
                    step = step + window
                    # if lyalphaind-step > -int(pw/2):
                    #     winpeakind = np.arange(step,(step+window))
                    #     winpeak = np.append(winpeak,flux[winpeakind])
                    #     intervalwlen = np.append(intervalwlen,wlen[winpeakind])
                    #     step = step + window

        #pad interval with start/end value to allign correctly
            #winpeakmed = step + findmed(windata)
            #startind = findmed(flux[0:window])
            winpeak[0] = flux[0]
            intervalwlen[0] = wlen[0]
            #endind = findmed(flux[-window:])
            winpeak[-1] = flux[-1]
            intervalwlen[-1] = wlen[-1]

            intpol = interpolate.interp1d(intervalwlen, winpeak, kind=1)
            contfit = intpol(wlen)
            normspec = flux-contfit


        #plotting:
            if len(specsample) < 10 and showplot == True:
                wlim = speclen
                plt.figure(spec[:20])
                #plt.title('continuum fitting')
                plt.plot(wlen[0:wlim],flux[0:wlim],label=spec[:20] + 'z = ' + str(redshift))
                plt.plot(intervalwlen[0:wlim],winpeak[0:wlim],'*',label='intervals')
                plt.plot(wlen[0:wlim],contfit[0:wlim],'--',label='interpolation')
                plt.plot(wlen[lyalphaind],flux[lyalphaind],'.',label=r'Ly$\alpha$')
                plt.xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
                plt.ylabel(r'$F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')
                #if specind == specsample[0]:
                plt.legend()

                plt.figure('continuum removed')
                #plt.title('continuum removed')
                plt.plot(wlen[0:wlim],normspec[0:wlim],label=spec[:20])
                plt.xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
                plt.ylabel(r'$F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')
                plt.legend()

    if showplot == True:
        plt.show()

    return wlen, normspec, lyalphacalc



#
#
# specsample = ['spec-0343-51692-0145.fits']
# for spec in specsample:
#     spec = [spec]
#     wlen, normspec, lyalpha = contfitv7(spec, showplot = True)
