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
#
# gcconf = ['J0116-2052' , 'J0800+4029','J0958−2904' ,'J1017+6116' ,'J1018+0530','J1052+0806',
# 'J1103+3449' ,'J1129+0951' ,'J1131−2705' ,'J1300+4009' ,'J1358+5752' ,'J1510+5958' ,
# 'J1753+6310' ,'J2039−2514' ,'J2227−2705' ,'J2355−0002']
#
# gcconfmatch = []
# for i in range(0,len(gcconf)):
#     sample = str(gcconf[i])
#     for testmatch in carlamatch:
#         trim = testmatch[0:5]
#         if trim == sample[0:5]:
#             gcconfmatch.append(testmatch)
# print(gcconfmatch)
# carlamatch = gcconfmatch #select just this subsample

########################################################


#do multiple bin stacks:
rinterval = 4000
rbins = np.array([0, rinterval])
while rbins[1] <= 4000:
    binrun = 'rad_bins_(' + str(rbins[0]) + '_to_' + str(rbins[1]) + ')'
    print('radial binning for ' + str(rbins[0]) + ' to ' + str(rbins[1]) + ' : ')

    specstacktot = 0 #total number of spectra stacked in all carla
    carlamatchlen = len(carlamatch)
    carlamatchlen = 10
    carlarange = np.arange(0,carlamatchlen)
    #variables for stacking carla together
    carlahighreslen = 100000
    wlenmultistack = np.linspace(500, 4000, carlahighreslen)
    carlacutinds = []
    meanmultistore = np.zeros([len(carlarange),carlahighreslen])
    medmultistore = np.zeros([len(carlarange),carlahighreslen])

    meancontmultistore = np.zeros([len(carlarange),carlahighreslen])
    medcontmultistore = np.zeros([len(carlarange),carlahighreslen])

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
        fillval = np.nan
        if len(specselect) == 0:#catched cluster
            print('')
            print('Cluster sample error!: '+ carlamatch[carlaselect] + ' has no spec that meet criteria' )
            print('')
            meanmultistore[carlanumber, 0:] = fillval
            medmultistore[carlanumber, 0:] = fillval
            meancontmultistore[carlanumber, 0:] = fillval
            medcontmultistore[carlanumber, 0:] = fillval
            carlacutinds.append(carlanumber)
        else:
            highreslen = 50000
            cutinds = []
            normspecstore = np.zeros([len(specselect),highreslen])
            contspecstore = np.zeros([len(specselect),highreslen])
            wlenhighres = np.linspace(500, 4500, highreslen)
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
                #get spec/fit data
                wlen = data.field(0)
                flux = data.field(1)
                medcontfit = data.field(2)
                maxcontfit = data.field(3)
                #continuum normalisation
                contfit = maxcontfit
                normspec = flux

                #get ivar from original spec
                oldspecdirectory = 'Spectra/'+spec[0:20]+'.fits'
                olddata = fits.getdata(oldspecdirectory,ext=1)
                var = 1/(olddata.field(3))

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
                zlims = np.array([gcredshift+0.05 , gcredshift + 4.5])
                stonlim = 1
                pw = 0
                #check fit can be added to stack:
                if gc_qso_sep < rbins[0] or gc_qso_sep > rbins[1]:
                    stackstatus.append('filtererror')
                    normspecstore[specnumber, 0:] = np.nan
                    contspecstore[specnumber, 0:] = np.nan
                    cutinds.append(specnumber)
                    if showerror == True:
                        print('filter warning!: '+spec+' z = '+ str(redshift) +' did not meet filter condition (S/N = '+ str(stonall) +')')
                elif lyalpha - wlen[0] <= 0: #ignore qso if it has no forest (lyalpha - wlen[0] <= 0)
                    stackstatus.append('foresterror')
                    normspecstore[specnumber, 0:] = np.nan
                    contspecstore[specnumber, 0:] = np.nan
                    cutinds.append(specnumber)
                    if showerror == True:
                        print('forest warning!: '+spec+' z = '+ str(redshift) +' forest before start of spectrum (S/N = '+ str(stonall) +')')
                elif redshift < zlims[0] or redshift > zlims[1]: #removes qso in carla (those within z = 0.05) and too high
                    stackstatus.append('zerror')
                    normspecstore[specnumber, 0:] = np.nan
                    contspecstore[specnumber, 0:] = np.nan
                    cutinds.append(specnumber)
                    if showerror == True:
                        print('zlim warning!: '+spec+' z = '+ str(redshift) +' below not in limits (S/N = '+ str(stonall) +')')
                elif stonall < stonlim:
                    stackstatus.append('stonerror')
                    normspecstore[specnumber, 0:] = np.nan
                    contspecstore[specnumber, 0:] = np.nan
                    cutinds.append(specnumber)
                    if showerror == True:
                        print('S/N warning!: '+spec+' S/N = '+ str(stonall) +' too low (z = '+ str(redshift) +')')
                elif gclyalpha - wlen[0] <= pw:
                    stackstatus.append('gcerror')
                    normspecstore[specnumber, 0:] = np.nan
                    contspecstore[specnumber, 0:] = np.nan
                    cutinds.append(specnumber)
                    if showerror == True:
                        print('gc warning!: '+spec+' has z too low with respect to gc withlyalpha '+ str(gclyalpha - wlen[0]) + ' Angstroms before start of spectra')
                else:
                    #continuum regulation as per Faucher and Giguierre
                    wlenshift = wlen/(1+gcredshift)

                    stackstatus.append('success')
                    #wlenshift = wlen/(1+redshift)
                    wlenintpol = interpolate.interp1d(wlenshift, normspec, 'linear', bounds_error=False, fill_value=fillval)
                    contintpol = interpolate.interp1d(wlenshift, contfit, 'linear', bounds_error=False, fill_value=fillval)
                    if wlenshift[0] < wlenmin:
                        wlenmin = wlenshift[0]
                    if wlenshift[-1] > wlenmax:
                        wlenmax = wlenshift[-1]

                    normspecstore[specnumber, 0:] = wlenintpol(wlenhighres)
                    contspecstore[specnumber, 0:] = contintpol(wlenhighres)

                    if showerror == True:
                        print('adding '+spec+' S/N = '+ str(stonall) +'and z = '+ str(redshift) +' to stack.')

                        #continuum plot check
                        range = np.arange(0,2500).astype(int)
                        fig1, ax = plt.subplots(2,1,num=spec[0:20])
                        ax[0].plot(wlenshift[range], flux[range], label='flux')
                        ax[0].plot(wlenshift[range], medcontfit[range], label='med')
                        ax[0].plot(wlenshift[range], maxcontfit[range], label='max')
                        ax[0].set_ylabel(r'$F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')
                        ax[0].legend()
                        ax[1].plot(wlenshift[range], normspec[range], label='normmax')
                        ax[1].plot(wlenshift[range], var[range], label='ivar')
                        ax[1].set_xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
                        gcabs = 1215.67
                        qsoabs = 1215.67/(1+gcredshift)*(1+redshift)
                        ax[1].plot(np.array([gcabs,gcabs]),np.array([0,1.5]),label='GC lyalpha abs')
                        ax[1].plot(np.array([qsoabs,qsoabs]),np.array([0,1.5]),label='QSO lyalpha em')
                        ax[1].legend()
                        #ax[1].set_ylim((-2, 2))   # set the ylim to bottom, top
                        plt.show()

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
                meanmultistore[carlanumber, 0:] = np.nan
                medmultistore[carlanumber, 0:] = np.nan
                contmultistore[carlanumber, 0:] = np.nan
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
                start = (np.abs(wlenhighres - wlenmin)).argmin() +10
                end = (np.abs(wlenhighres - wlenmax)).argmin() -10
                wlenhighres = wlenhighres[start:end]
                #cut downcols / remove empty rows
                #normspecstore = np.delete(normspecstore, cutinds, axis = 0)
                normspecstore = normspecstore[0:, start:end]
                meanspec = np.nanmean(normspecstore, axis=0)
                medspec = np.nanmedian(normspecstore, axis=0)

                contspecstore = contspecstore[0:, start:end]
                meancont = np.nanmean(contspecstore, axis=0)
                medcont = np.nanmedian(contspecstore, axis=0)


                #stack carla
                meangcwlenintpol = interpolate.interp1d(wlenhighres, meanspec, 'linear', bounds_error=False, fill_value=fillval)
                medgcwlenintpol = interpolate.interp1d(wlenhighres, medspec, 'linear', bounds_error=False, fill_value=fillval)
                meangccontintpol = interpolate.interp1d(wlenhighres, meancont, 'linear', bounds_error=False, fill_value=fillval)
                medgccontintpol = interpolate.interp1d(wlenhighres, medcont, 'linear', bounds_error=False, fill_value=fillval)
                if wlenhighres[0] < gcwlenmin:
                    gcwlenmin = wlenhighres[0]
                if wlenhighres[-1] > gcwlenmax:
                    gcwlenmax = wlenhighres[-1]

                #append stack to multidim variable:
                meanmultistore[carlanumber, 0:] = meangcwlenintpol(wlenmultistack)
                medmultistore[carlanumber, 0:] = medgcwlenintpol(wlenmultistack)

                meancontmultistore[carlanumber, 0:] = meangccontintpol(wlenmultistack)
                medcontmultistore[carlanumber, 0:] = medgccontintpol(wlenmultistack)

        carlanumber = carlanumber + 1
    #print uncombined stacks for carla


    #cut down zero padding
    carlasucces = carlanumber - len(carlacutinds)
    gcstart = (np.abs(wlenmultistack - gcwlenmin)).argmin() +10 #tolerance for error in exact start
    gcend = (np.abs(wlenmultistack - gcwlenmax)).argmin() -10
    wlenmultistack = wlenmultistack[gcstart:gcend]

    #cut downcols / remove empty rows
    #cut carlamatch

    #mean
    meanmultistore = meanmultistore[0:, gcstart:gcend]
    meancarla = np.nanmean(meanmultistore, axis=0)

    meancontmultistore = meancontmultistore[0:, gcstart:gcend]
    meancontcarla = np.nanmean(meancontmultistore, axis=0)

    #med
    medmultistore = medmultistore[0:, gcstart:gcend]
    medcarla = np.nanmedian(medmultistore, axis=0)

    medcontmultistore = medcontmultistore[0:, gcstart:gcend]
    medcontcarla = np.nanmean(medcontmultistore, axis=0)

    plt.figure('raw mean')
    plt.plot(wlenmultistack,meancarla)
    plt.plot(wlenmultistack,meancontcarla)



    #continuum regulation as per Faucher and Giguierre
    range = np.arange(findval(wlenmultistack, 1100),findval(wlenmultistack, 3000))
    wlen = wlenmultistack[range]*(1+gcredshift)
    wlenshift = wlenmultistack[range]
    flux = meancarla[range]
    cont = meancontcarla[range]

    norm = flux/cont

    lya = 1215.67
    z = (wlen/lya) -1
    zshift = z-gcredshift
    a = 0.0018
    b = 3.92
    teff = a*(1+z)**b
    faucher = np.exp(-teff)

    #fit norm continuum
    def quad(x, a, b, c):
        return a*x**2+b*x+c

    popt, pcov = cf(quad, wlenshift, norm)
    a, b, c = popt
    fit = quad(wlenshift,a, b, c)
    corrected =  norm + (faucher - fit)


    fig1, ax = plt.subplots(2,1,num='plot')

    ax[0].plot(wlenshift, flux, label='flux')
    ax[0].plot(wlenshift, cont, label='cont')
    ax[0].set_xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
    ax[0].set_ylabel(r'$<F>$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')
    ax[0].legend()

    ax[1].plot(wlenshift, (norm), label='norm')
    ax[1].plot(wlenshift, (faucher), label='faucher')
    ax[1].plot(wlenshift, fit, label='fit')
    ax[1].plot(wlenshift, corrected, label='corrected')
    ax[1].set_xlabel(r'$z$')
    ax[1].legend()
    plt.show()




    rbins[0] = rbins[0] + rinterval
    rbins[1] = rbins[1] + rinterval
