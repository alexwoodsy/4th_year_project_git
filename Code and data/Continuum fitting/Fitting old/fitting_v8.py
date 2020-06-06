#top x% method for fitting continuum in the forest
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
import os

#plt.style.use('mystyle') #path C:\Users\alexw\AppData\Local\Programs\Python\Python37\Lib\site-packages\matplotlib\mpl-data\stylelib

#calculates points slected in interval
def findval(array,val):
    array = np.asarray(array)
    ind = np.abs(array - val).argmin()
    return ind

def findmax(array):
    array = np.asarray(array)
    ind = (np.abs(array - np.max(array))).argmin()
    return ind

def findmed(array):
    array = np.asarray(array)
    ind = (np.abs(array - np.median(array))).argmin()
    return ind

def findpctmax(array,pct):
    ind = np.argsort(array)
    sortedarray = array[ind]
    selectedinds = ind[-pct:] #indices of the top pct in array
    selectedinds = ind[int(-pct*0.5)]
    return selectedinds

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

#ref lines from sdss
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


def contfitv8(specsample,zlim,stonlim,gcredshift,showplot,showerror):

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



    #z correct refline and find corresponding index 9 (finds ly alpha too)
        wlenlinezcorr = wlenline*(1+redshift)
        wlenlineind = np.array([]).astype(int)
        for i in range(0,len(wlenlinezcorr)):
            ind = (np.abs(wlen - wlenlinezcorr[i])).argmin()#finds index of nearest point in data
            wlenlineind = np.append(wlenlineind,ind).astype(int)

    #lyalpha considerations
        lyalpha = 1215.67*(1+redshift)#calc lya using redshift
        gclyalpha = 1215.67*(1+gcredshift)
        lyalphaind = wlenlineind[1]


    #fitting:
        #split the spec in two about lyalpha peak
        peaklim = 0 #zero to fit peaks
        pw = 40# for removinf ly alpha peak
        intervalwlen = np.array([])
        winpeak = np.array([])
        forestwinnum, forestpct = 20, 0.4 #forest number must be even
        otherwinnum, otherpct =  50, 0.2


        #s/n check and
        #loop increments print('z warning!: '+spec+' has z too low with  lyalpha '+str(lyalpha - wlen[0]+' Anstroms before start of spectra')

        if redshift < zlim and showerror:
            normspec = np.zeros(speclen)
            stackstatus = 'zerror'
            if showerror == True:
                print('z warning!: '+spec+' z = '+ str(redshift) +' too low (S/N = '+ str(ston) +') - normspec = 0 array')
        elif ston < stonlim:
            normspec = np.zeros(speclen)
            stackstatus = 'stonerror'
            if showerror == True:
                print('S/N warning!: '+spec+' S/N = '+ str(ston) +' too low (z = '+ str(redshift) +') - normspec = 0 array')
        elif gclyalpha - wlen[0] <= pw:
            normspec = np.zeros(speclen)
            stackstatus = 'foresterror'
            if showerror == True:
                print('z warning!: '+spec+' has z too low with respect to gc withlyalpha '+ str(lyalpha - wlen[0]) + ' Angstroms before start of spectra')
        else:
            #proceed with continuum fitting
            stackstatus = 'success'
            if showerror == True:
                print('adding '+spec+' S/N = '+ str(ston) +'and z = '+ str(redshift) +' to stack.')
            step = 0
            while step <= speclen:
                if step <= lyalphaind:
                    winnum = forestwinnum
                    window = int(lyalphaind/winnum)
                    percentage = forestpct
                    pct = int(window*percentage)
                    windata = flux[step:(step+window)]
                    winpeakind = step + findpctmax(windata,pct)#change
                    winpeak = np.append(winpeak,flux[winpeakind])
                    intervalwlen = np.append(intervalwlen,wlen[winpeakind])

                    for ind in range(step,(step+window)): #removes ly alpha
                        if lyalphaind-ind == pw: #ensures no dip before remving lyalpha
                            startind = step - int(pw/2) + findmax(flux[(step-int(pw/2)):ind])
                            winpeak = np.append(winpeak,flux[startind])
                            intervalwlen = np.append(intervalwlen,wlen[startind])
                        if abs(lyalphaind-ind < pw):
                            winpeak = np.append(winpeak,flux[ind])
                            intervalwlen = np.append(intervalwlen,wlen[ind])


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
            smoothwin = 4*int(lyalphaind/forestwinnum)+1 # has to be odd hence forestwinnum must be even
            #print(smoothwin)
            forrange = np.arange(0,(lyalphaind-pw-1))
            contfit[forrange] = signal.savgol_filter(contfit[forrange], smoothwin,3)
            normspec = flux-contfit





        #plotting:
            if len(specsample) < 10 and showplot == True:
                wlim = speclen
                plt.figure(spec[:20])
                #plt.title('continuum fitting')
                plt.plot(wlen[0:wlim],flux[0:wlim],label=spec[:20] + 'z = ' + str(redshift))
                plt.plot(intervalwlen[0:wlim],winpeak[0:wlim],'*',label='intervals')
                plt.plot(wlen[0:wlim],contfit[0:wlim],'--',label='interpolation')
                plt.plot(wlen[wlenlineind],flux[wlenlineind],'.',label='linerefs')
                for i in range(0,len(wlenlineind)):
                    plt.text(wlen[wlenlineind[i]], flux[wlenlineind[i]]+0.3, linename[i], fontsize=9)
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

    #total errors and total stacked


    return wlen, normspec, wlenlineind, redshift, stackstatus



#
#
# specsample = ['spec-0343-51692-0145.fits']
# for spec in specsample:
#     spec = [spec]
#     wlen, normspec, lyalpha = contfitv7(spec, showplot = True)
