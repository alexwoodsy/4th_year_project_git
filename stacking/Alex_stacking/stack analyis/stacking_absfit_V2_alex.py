#stacking v2
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
from scipy.optimize import curve_fit as cf
import os, time
#import fitting
import sys

#plt.style.use('mystyle') #path C:\Users\alexw\AppData\Local\Programs\Python\Python37\Lib\site-packages\matplotlib\mpl-data\stylelib

def findval(array,val):
    array = np.asarray(array)
    ind = np.abs(array - val).argmin()
    return ind

#imports the  FITTED spectra from the spectra folder
specnames = next(os.walk('Fitted Spectra'))[2]
spectot = len(specnames)

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

#get number of carla targets in sample
#read in carla
carladata = fits.getdata('CARLA/CARLA_table_Aug_2013.fits',ext=1)
carlalen = len(carladata)
carlanames = np.zeros(carlalen).astype(str)
overdensity = np.zeros(carlalen)
carlamatch =  []
carlaconfirmed = 1
for i in range(0,carlalen):
    c = 0
    carlanames = str(carladata[i][0])
    overdensity[i] = carladata[i][6]
    for k in range(0,fitlen):
        if carlanames ==  metagcname[k]: #'J080016.10+402955.6':#and overdensity[i] > 3: #add conditional for overdendity etc.. here
            c = c + 1
    if c > 0:
        carlamatch.append(carlanames)

##### TEST SAMPLE OF CONFIRMED GCS AROUND CARLA targets

gcconf = ['J0116-2052' , 'J0800+4029','J0958−2904' ,'J1017+6116' ,'J1018+0530','J1052+0806',
'J1103+3449' ,'J1129+0951' ,'J1131−2705' ,'J1300+4009' ,'J1358+5752' ,'J1510+5958' ,
'J1753+6310' ,'J2039−2514' ,'J2227−2705' ,'J2355−0002']

gcconfmatch = []
for i in range(0,len(gcconf)):
    sample = str(gcconf[i])
    for testmatch in carlamatch:
        trim = testmatch[0:5]
        if trim == sample[0:5]:
            gcconfmatch.append(testmatch)
print(gcconfmatch)
carlamatch = gcconfmatch #select just this subsample

########################################################


#do multiple bin stacks:
rinterval = 4000
rbins = np.array([0, rinterval])
while rbins[1] <= 4000:
    binrun = 'rad_bins_(' + str(rbins[0]) + '_to_' + str(rbins[1]) + ')'
    print('radial binning for ' + str(rbins[0]) + ' to ' + str(rbins[1]) + ' : ')

    specstacktot = 0 #total number of spectra stacked in all carla
    carlamatchlen = len(carlamatch)
    carlarange = np.arange(0,carlamatchlen)
    #variables for stacking carla together
    carlahighreslen = 100000
    wlenmultistack = np.linspace(500, 4000, carlahighreslen)
    carlacutinds = []
    meanmultistore = np.zeros([len(carlarange),carlahighreslen])
    medmultistore = np.zeros([len(carlarange),carlahighreslen])
    gcwlenmin = carlahighreslen
    gcwlenmax = 0

    carlanumber = 0
    for carlaselect in carlarange: #all change to - matchlen
        print('processing '+ carlamatch[carlaselect]+ ' ('+str(carlaselect)+'/'+str(len(carlarange))+')')

        #get sample to stack:
        specselect = []
        for k in range(0,fitlen):
            if metagcname[k] == carlamatch[carlaselect]:  #CHNAGE DEPEDNING ON SAMPLE
                specselect.append(k)

        print('sample of '+str(len(specselect))+' found for stacking')
        print('|')

        #####----STACKING---#####
        fillval = 1
        if len(specselect) == 0:#catched cluster
            print('')
            print('Cluster sample error!: '+ carlamatch[carlaselect] + ' has no spec that meet criteria' )
            print('')
            meanmultistore[carlanumber, 0:] = np.zeros(carlahighreslen)
            medmultistore[carlanumber, 0:] = np.zeros(carlahighreslen)
            carlacutinds.append(carlanumber)
        else:
            highreslen = 50000
            cutinds = []
            normspecstore = np.zeros([len(specselect),highreslen])
            wlenhighres = np.linspace(500, 4000, highreslen)
            wlenmin = 10000
            wlenmax = 0
            stackstatus = []

            #debugging parameters:
            showerror = False
            specnumber = 0

            for ind in specselect:
                spec = metaspecfilename[ind]
                #import spectrum data:
                specdirectory = 'Fitted Spectra/'+spec
                data = fits.getdata(specdirectory,ext=1)
                speclen = len(data)

                wlen = data.field(0)
                flux = data.field(1)
                medcontfit = data.field(2)
                maxcontfit = data.field(3)

                normspec = flux/medcontfit

                #extract metadata
                redshift= metaredshift[ind]
                stonall = metastonall[ind]
                stonforest = metastonforest[ind]
                lyalpha = metalyalpha[ind]
                lyalphaind = metalyalphaind[ind]
                gcname = metagcname[ind]
                gcredshift = metagcredshift[ind]
                gc_qso_sep = metagc_qso_sep[ind]
                gclyalpha = metagclyalpha[ind]
                gclyalphaind = metagclyalphaind[ind]
                stackmsg = metastackmsg[ind]


                #conditions for stacking
                zlims = np.array([gcredshift+0.1 , gcredshift + 5])
                stonlim = 1
                pw = 0
                #check fit can be added to stack:
                if gc_qso_sep < rbins[0] or gc_qso_sep > rbins[1]:
                    stackstatus.append('filtererror')
                    normspecstore[specnumber, 0:] = np.zeros(highreslen)
                    cutinds.append(specnumber)
                    if showerror == True:
                        print('filter warning!: '+spec+' z = '+ str(redshift) +' did not meet filter condition (S/N = '+ str(stonall) +')')
                elif lyalpha - wlen[0] <= 0: #ignore qso if it has no forest (lyalpha - wlen[0] <= 0)
                    stackstatus.append('foresterror')
                    normspecstore[specnumber, 0:] = np.zeros(highreslen)
                    cutinds.append(specnumber)
                    if showerror == True:
                        print('forest warning!: '+spec+' z = '+ str(redshift) +' forest before start of spectrum (S/N = '+ str(stonall) +')')
                elif redshift < zlims[0] or redshift > zlims[1]: #removes qso in carla (those within z = 0.05) and too high
                    stackstatus.append('zerror')
                    normspecstore[specnumber, 0:] = np.zeros(highreslen)
                    cutinds.append(specnumber)
                    if showerror == True:
                        print('zlim warning!: '+spec+' z = '+ str(redshift) +' below not in limits (S/N = '+ str(stonall) +')')
                elif stonall < stonlim:
                    stackstatus.append('stonerror')
                    normspecstore[specnumber, 0:] = np.zeros(highreslen)
                    cutinds.append(specnumber)
                    if showerror == True:
                        print('S/N warning!: '+spec+' S/N = '+ str(stonall) +' too low (z = '+ str(redshift) +')')
                elif gclyalpha - wlen[0] <= pw:
                    stackstatus.append('gcerror')
                    normspecstore[specnumber, 0:] = np.zeros(highreslen)
                    cutinds.append(specnumber)
                    if showerror == True:
                        print('gc warning!: '+spec+' has z too low with respect to gc withlyalpha '+ str(gclyalpha - wlen[0]) + ' Angstroms before start of spectra')
                else:
                    #proceed with continuum fitting
                    stackstatus.append('success')
                    if showerror == True:
                        print('adding '+spec+' S/N = '+ str(stonall) +'and z = '+ str(redshift) +' to stack.')

                    wlenshift = wlen/(1+gcredshift)
                    wlenintpol = interpolate.interp1d(wlenshift, normspec, 'linear', bounds_error=False, fill_value=fillval)
                    if wlenshift[0] < wlenmin:
                        wlenmin = wlenshift[0]
                    if wlenshift[-1] > wlenmax:
                        wlenmax = wlenshift[-1]
                    normspechighres = wlenintpol(wlenhighres)
                    normspecstore[specnumber, 0:] = normspechighres
                specnumber = specnumber + 1
            #output
            print('stacking attempted for '+ str(len(stackstatus)) + ' spectra for CARLA: ' + carlamatch[carlaselect])

            stacktot = zerrortot = foresterrortot = stonerrortot = gcerrortot =filtererrortot = 0

            for x in stackstatus:
                if x == 'success':
                    stacktot = stacktot + 1
                if x == 'gcerror':
                    gcerrortot = gcerrortot + 1
                if x == 'zerror':
                    zerrortot = zerrortot + 1
                if x == 'foresterror':
                    foresterrortot = foresterrortot + 1
                if x =='stonerror':
                    stonerrortot = stonerrortot + 1
                if x =='filtererror':
                    filtererrortot = filtererrortot + 1

            errortot = zerrortot + foresterrortot + stonerrortot + gcerrortot + filtererrortot

            if errortot == len(stackstatus):
                print('stack warning!: no spec stacked, with ' + str(errortot) + ' not used:')
                print('')
                meanmultistore[carlanumber, 0:] = np.zeros(carlahighreslen)
                medmultistore[carlanumber, 0:] = np.zeros(carlahighreslen)
                carlacutinds.append(carlanumber)
            else:
                print(str(stacktot) +' stacked successfully with ' + str(errortot) + ' not used:')
                print(str(filtererrortot) + ' didn not meet filter conditon ')
                print(str(zerrortot) + ' outside redshift range ' + str(zlims[0]) + ' < z < '+ str(zlims[1]))
                print(str(foresterrortot) + ' lyalpaforest out of spectra range' )
                print(str(gcerrortot) + ' galaxy luxter lyalpha out of spectra range' )
                print(str(stonerrortot) + ' S/N below ' + str(stonlim))
                print(' ')

                specstacktot = specstacktot + stacktot

                #stacking data only if there are spec to stack
                #cut extra zeropadding
                start = (np.abs(wlenhighres - wlenmin)).argmin()
                end = (np.abs(wlenhighres - wlenmax)).argmin()
                wlenhighres = wlenhighres[start:end]
                #cut downcols / remove empty rows
                normspecstore = np.delete(normspecstore, cutinds, axis = 0)
                normspecstore = normspecstore[0:, start:end]
                meanspec = np.mean(normspecstore, axis=0)
                medspec = np.median(normspecstore, axis=0)


                #stack carla
                meangcwlenintpol = interpolate.interp1d(wlenhighres, meanspec, 'linear', bounds_error=False, fill_value=fillval)
                medgcwlenintpol = interpolate.interp1d(wlenhighres, medspec, 'linear', bounds_error=False, fill_value=fillval)
                if wlenhighres[0] < gcwlenmin:
                    gcwlenmin = wlenhighres[0]
                if wlenhighres[-1] > gcwlenmax:
                    gcwlenmax = wlenhighres[-1]
                meancarlastack = meangcwlenintpol(wlenmultistack)
                medcarlastack = medgcwlenintpol(wlenmultistack)
                #append stack to multidim variable:
                meanmultistore[carlanumber, 0:] = meancarlastack
                medmultistore[carlanumber, 0:] = medcarlastack

        carlanumber = carlanumber + 1
    #print uncombined stacks for carla


    #cut down zero padding
    carlasucces = carlanumber - len(carlacutinds)
    gcstart = (np.abs(wlenmultistack - gcwlenmin)).argmin()
    gcend = (np.abs(wlenmultistack - gcwlenmax)).argmin()
    wlenmultistack = wlenmultistack[gcstart:gcend]

    #cut downcols / remove empty rows
    #cut carlamatch

    #mean
    meanmultistore = np.delete(meanmultistore, carlacutinds, axis = 0)
    meanmultistore = meanmultistore[0:, gcstart:gcend]
    meancarla = np.mean(meanmultistore, axis=0)
    #med
    medmultistore = np.delete(medmultistore, carlacutinds, axis = 0)
    medmultistore = medmultistore[0:, gcstart:gcend]
    medcarla = np.median(medmultistore, axis=0)


    ###absorbtion line fitting####
    #fitting needs initial data so extract data about abs line
    def guassian(x, amp, mean, std):
        return amp*np.exp(-((x-mean)**2)/(2*(std)**2))

    fig0, ax = plt.subplots(2,1,num=binrun+'Absorption line fitting')

    abslineind = findval(wlenmultistack, 1215.67)
    datarange = np.arange(abslineind - 12, abslineind + 15)
    plotrange = np.arange(abslineind - 50, abslineind + 50)
    absflux = meancarla[datarange]
    abswlen = wlenmultistack[datarange]

    popt, pcov = cf(guassian, abswlen, absflux, bounds =([-np.inf,1215.67-0.5,-np.inf],[np.inf,1215.67+0.5,np.inf]))
    meanamp, meanmean, meanstd = popt
    ax[0].plot(wlenmultistack[plotrange], meancarla[plotrange])
    ax[0].plot(abswlen, guassian(abswlen, *popt), 'r-',label='fitting parmaters: amp=%5.3f, mean=%5.3f, std=%5.3f' % tuple(popt))
    ax[0].set_xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
    ax[0].set_ylabel(r'$<F>$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')
    ax[0].legend()

    absflux = medcarla[datarange]
    abswlen = wlenmultistack[datarange]
    popt, pcov = cf(guassian, abswlen, absflux, bounds =([-np.inf,1215.67-0.5,-np.inf],[np.inf,1215.67+0.5,np.inf]))
    medamp,medmean,medstd = popt
    ax[1].plot(wlenmultistack[plotrange], medcarla[plotrange])
    ax[1].plot(abswlen, guassian(abswlen, *popt), 'r-',label='fitting parmaters: amp=%5.3f, mean=%5.3f, std=%5.3f' % tuple(popt))
    ax[1].set_xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
    ax[1].set_ylabel(r'MEDIAN $F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')
    ax[1].legend()

    #plot relative vel graph and find delv of line from rest gc ly alpha
    fig1, ax = plt.subplots(2,1,num=binrun+'velocity Absorption line plot')
    c = 299792.458
    lam = wlenmultistack
    lam_em = 1215.67
    vrel = c*((lam  - lam_em)/lam_em)

    abslineind = findval(vrel, 0)
    plotrange = np.arange(abslineind - 50, abslineind + 50)

    ax[0].plot(vrel, meancarla)
    ax[0].set_xlabel(r'$\delta$v ($kms^{-1})$')
    ax[0].set_ylabel(r'$<F>$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')

    ax[1].plot(vrel, medcarla)
    ax[1].set_xlabel(r'$\delta$v ($kms^{-1}$)')
    ax[1].set_ylabel(r'MEDIAN $F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')


    # plt.figure() #voigt profile
    # x = np.arange(0, 10, 0.01)
    # v1 = models.Voigt1D(x_0=5, amplitude_L=10, fwhm_L=0.5, fwhm_G=0.9)
    # plt.plot(x, v1(x))

    #plotting output
    fig2, ax = plt.subplots(2,1,num=binrun+'multi-carla stack - uncombined')
    for c in range(0,carlasucces):

        dsrange = np.linspace(wlenmultistack[0], wlenmultistack[-300],5000)
        dsflux = signal.resample(meanmultistore[c, :-300], 5000)
        ax[0].plot(dsrange, dsflux)
        ax[0].set_xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
        ax[0].set_ylabel(r'$<F>$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')

        dsflux = signal.resample(medmultistore[c, :-300], 5000)
        ax[1].plot(dsrange, dsflux)
        ax[1].set_xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
        ax[1].set_ylabel(r'MEDIAN $F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')

    fig3, ax = plt.subplots(2,1,num=binrun+'('+str(carlasucces)+'/'+str(carlanumber)+ ') stacked clusters')
    ax[0].plot(wlenmultistack[:-300], meancarla[:-300])
    ax[0].plot(np.array([(1215.67),(1215.67)]),np.array([np.min(meancarla),np.max(meancarla)]),'--')
    ax[0].text(1215.67, np.min(meancarla), r' Ly$\alpha$ absorbtion for multi-carla stack with '+ str(specstacktot) +' spectra')
    ax[0].set_xlabel(r' $\lambda$ ($\mathrm{\AA}$)')
    ax[0].set_ylabel(r'$<F>$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')

    ax[1].plot(wlenmultistack[:-300], medcarla[:-300])
    ax[1].plot(np.array([(1215.67),(1215.67)]),np.array([np.min(medcarla),np.max(medcarla)]),'--')
    ax[1].text(1215.67, np.min(medcarla), r' Ly$\alpha$ absorbtion for multi-carla stack with '+ str(specstacktot) +' spectra')
    ax[1].set_xlabel(r' $\lambda$ ($\mathrm{\AA}$)')
    ax[1].set_ylabel(r'MEDIAN $F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')

    plt.show()

    rbins[0] = rbins[0] + rinterval
    rbins[1] = rbins[1] + rinterval
