#---task overview
#loop through all spec in spectra folder
#filter spectra accordingly so that those we cannot fit are ignored (maybe pass through fit in bit we cant fit)
#perform fitting by fitting_v9 + slighlty modified fitting with alternate window sizes
#nt - fit algorthm must be fully optimised e.g, removed figures / unneccsary steps
#---jason save code
#----metadata to save:
# s/n calculated for whole spec + forest
#associated galaxy cluster + its gcredshift
#redshift of qso
#position of gclyalpha in spec (-ve means before start of spec! )
#seperation from gc in arcsecs on the sky
#fitting data = window size of savgol_filter and order
#flag to state whether there was a forest to even fit properly

import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
import os


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


#imports the spectra from the spectra folder
specnames = next(os.walk('Spectra'))[2]
spectot = len(specnames)
#add indexing for spectra in file to allow loop over all
#specnames = [specnames[0],specnames[1000]]


for spec in specnames:

    #---------------------------data extraction-----------------------------#
    specdirectory = 'Spectra/'+spec
    #print(specdirectory)

    data = fits.getdata(specdirectory,ext=1)#Read in fits spectrum data
    speclen = len(data)
    fitdata = fits.getdata(specdirectory,ext=2)#read in fits meta data
    metasize = len(fitdata[0])
    #predefine data variables
    flux = np.zeros(speclen)
    wlen = np.zeros(speclen)
    model = np.zeros(speclen)
    ivar = np.zeros(speclen)

    #extract metadata
    for i in range(0,speclen):
        flux[i] = data[i][0]
        ivar[i] = data[i][2]
        if ivar[i] == 0:
            ivar[i] = 0.0000001
        wlen[i] = 10**(data[i][1])
        model[i] = data[i][7]

    #extract redshift of qso
    if metasize == 126:
         redshift = fitdata[0][63] #accounts for diff location for sdss/eboss
    else:
         redshift = fitdata[0][38]

    #extract qso-gc info from positions table -> asscoaited gc and angular seperation etc...
    posdata = fits.getdata('PositionsTable.fits',ext=1)
    poslen = len(posdata)
    clusterredshift = 0 #metadata
    clusternames = 'null' #metadata
    clusterseperation = 0 #metadata

    for j in range(0,poslen):
        #get assocaited id to match spec in folder to spec in postions table
        plate = str(posdata[j][4]).zfill(4)
        mjd = str(posdata[j][5]).zfill(5)
        fiberid = str(posdata[j][6]).zfill(4)
        specfilename = 'spec-'+plate+'-'+mjd+'-'+fiberid+'.fits'
        if specfilename == spec:
            clusterredshift = posdata[j][109]
            clusternames = str(posdata[j][105])
            clusterseperation = posdata[j][114]


    #lymanalpha calculation:
    lyalpha = 1215.67*(1+redshift)#calc lya using redshift #metadata
    lyalphaind = findval(wlen,lyalpha) #metadata
    gclyalpha = 1215.67*(1+clusterredshift) #metadata
    gclyalphaind = findval(wlen,gclyalpha) #metdata

    #s/n checking
    std = (1/ivar)**0.5
    stonall = np.median(flux/std) #metadata
    stonforest = np.median(flux[0:lyalphaind]/std[0:lyalphaind]) #metadata

    #--------------------Continuum fitting-----------------------#

    #split the spec in two about lyalpha peak for 2 fitting regions
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
        normspec = flux-contfit

        wlencol = fits.Column(name='Wavelength', array=wlen, format='F')
        normspeccol = fits.Column(name='Normalised_Flux', array=normspec, format='F')

        specfitdata = fits.BinTableHDU.from_columns([wlencol, normspeccol])

        redshiftcol = fits.Column(name='QSO_Z', array=np.array([redshift]), format='F')
        stonallcol = fits.Column(name='STON_ALL', array=np.array([stonall]), format='F')
        stonforestcol = fits.Column(name='STON_FOREST', array=np.array([stonforest]), format='F')
        lyalphacol = fits.Column(name='QSO_LYa', array=np.array([lyalpha]), format='F')
        lyalphaindcol = fits.Column(name='QSO_LYa_INDEX', array=np.array([lyalphaind]), format='K')
        gcnamecol = fits.Column(name='GC_NAME', array=np.array([clusternames]), format='20A')
        gcredshiftcol = fits.Column(name='GC_Z', array=np.array([clusterredshift]), format='F')
        gc_qso_sepcol = fits.Column(name='GC_SEPERATION', array=np.array([clusterseperation]), format='F')
        gclyalphacol = fits.Column(name='GC_LYa', array=np.array([gclyalpha]), format='F')
        gclyalphaindcol = fits.Column(name='GC_LYa_INDEX', array=np.array([gclyalphaind]), format='K')
        stackmsgcol = fits.Column(name='STACK_STATUS', array=np.array([stackstatus]), format='20A')

        metadata = fits.BinTableHDU.from_columns([redshiftcol, stonallcol, stonforestcol, lyalphacol,
        lyalphaindcol, gcnamecol, gcredshiftcol, gc_qso_sepcol, gclyalphacol, gclyalphaindcol, stackmsgcol])

        primary = fits.PrimaryHDU()
        hdul = fits.HDUList([primary, specfitdata, metadata])

        hdul.writeto('Fitted Spectra/' + spec[0:20] + '-prefitted.fits',overwrite = True)

        #plotting:
        # wlim = speclen
        # plt.figure(spec[:20]+'fitting')
        # plt.plot(wlen[0:wlim],flux[0:wlim],label=spec[:20] + 'z = ' + str(redshift))
        # plt.plot(intervalwlen[0:wlim],winpeak[0:wlim],'*',label='intervals')
        # plt.plot(wlen[0:wlim],contfit[0:wlim],'--',label='interpolation smoothed')
        # plt.xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
        # plt.ylabel(r'$F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')
        # plt.legend()
        #
        # plt.figure('continuum removed')
        # #plt.title('continuum removed')
        # plt.plot(wlen[0:wlim],normspec[0:wlim],label=spec[:20])
        # plt.xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
        # plt.ylabel(r'$F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')
        # plt.legend()
        # plt.show()






















    #
