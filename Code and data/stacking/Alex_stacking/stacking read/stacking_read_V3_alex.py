#stacking v2
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
import os, time
#import fitting
import sys
sys.path.append('C:/Users/alexw/Documents/GitHub/4th_year_project_git/Continuum fitting')
#path for other pc sys.path.append('C:/Users/alexw/OneDrive/Documents/University work/4th year work/Main project/4th_year_project_git/Continuum fitting')
import fitting_v9 as fitmeth
plt.style.use('mystyle') #path C:\Users\alexw\AppData\Local\Programs\Python\Python37\Lib\site-packages\matplotlib\mpl-data\stylelib

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
carlamatch =  []
for i in range(0,carlalen):
    c = 0
    carlanames = str(carladata[i][0])
    for k in range(0,fitlen):
        if carlanames == metagcname[k]: #add conditional for overdendity etc.. here
            c = c + 1
    if c > 0:
        carlamatch.append(carlanames)

#variables for stacking carla together
wlenmultistack = np.linspace(500, 4000, 100000)
multistack = np.zeros(100000)
gcwlenmin = 10000
gcwlenmax = 0

specstacktot = 0 #total number of spectra stacked in all carla
carlatot = 0
carlarange = np.arange(0,100)

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


    if len(specselect) == 0:#catched cluster
        print('')
        print('Cluster sample error!: '+ carlamatch[carlaselect] + ' has no spec that meet criteria' )
        print('')

    else:
        normspeckstack = np.zeros(50000)
        wlenhighres = np.linspace(500, 4000, 50000)
        wlenmin = 10000
        wlenmax = 0
        stackstatus = []

        #debugging parameters:
        showerror = False
        for ind in specselect:
            spec = metaspecfilename[ind]
            #import spectrum data:
            specdirectory = 'Fitted Spectra/'+spec
            data = fits.getdata(specdirectory,ext=1)
            speclen = len(data)

            normspec = data.field(1)
            wlen = data.field(0)

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
            zlims = np.array([gcredshift + 0.2, gcredshift + 5])
            stonlim = 1
            pw = 0

            #check fit can be added to stack:
            if lyalpha - wlen[0] <= 0: #ignore qso if it has no forest (lyalpha - wlen[0] <= 0)
                stackstatus.append('foresterror')
                if showerror == True:
                    print('forest warning!: '+spec+' z = '+ str(redshift) +' forest before start of spectrum (S/N = '+ str(stonall) +')')
            elif redshift < zlims[0] or redshift > zlims[1]: #removes qso in carla (those within z = 0.05) and too high
                stackstatus.append('zerror')
                if showerror == True:
                    print('zlim warning!: '+spec+' z = '+ str(redshift) +' below not in limits (S/N = '+ str(stonall) +')')
            elif stonall < stonlim:
                stackstatus.append('stonerror')
                if showerror == True:
                    print('S/N warning!: '+spec+' S/N = '+ str(stonall) +' too low (z = '+ str(redshift) +')')
            elif gclyalpha - wlen[0] <= pw:
                stackstatus.append('gcerror')
                if showerror == True:
                    print('gc warning!: '+spec+' has z too low with respect to gc withlyalpha '+ str(gclyalpha - wlen[0]) + ' Angstroms before start of spectra')
            else:
                #proceed with continuum fitting
                stackstatus.append('success')
                if showerror == True:
                    print('adding '+spec+' S/N = '+ str(stonall) +'and z = '+ str(redshift) +' to stack.')


                wlenshift = wlen/(1+gcredshift)
                wlenintpol = interpolate.interp1d(wlenshift, normspec, 'linear', bounds_error=False, fill_value=0)
                if wlenshift[0] < wlenmin:
                    wlenmin = wlenshift[0]
                if wlenshift[-1] > wlenmax:
                    wlenmax = wlenshift[-1]
                normspechighres = wlenintpol(wlenhighres)
                normspeckstack = normspeckstack + normspechighres

        #output
        print('stacking attempted for '+ str(len(stackstatus)) + ' spectra for CARLA: ' + carlamatch[carlaselect])

        stacktot = zerrortot = foresterrortot = stonerrortot = gcerrortot = 0

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

        errortot = zerrortot + foresterrortot + stonerrortot + gcerrortot

        if errortot == len(stackstatus):
            print('stack warning!: no spec stacked, with ' + str(errortot) + ' not used:')
            print('')
        else:
            carlatot = carlatot + 1
            print(str(stacktot) +' stacked successfully with ' + str(errortot) + ' not used:')
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
            normspeckstack = normspeckstack[start:end]

            plt.figure('multi-carla stack + uncombined')
            plt.plot(wlenhighres[:-300], normspeckstack[:-300], label = carlamatch[carlaselect] + ' stack of '+str(stacktot) + ' QSO')
            plt.xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
            plt.ylabel(r'$F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')

            #stack carla
            gcwlenintpol = interpolate.interp1d(wlenhighres, normspeckstack, 'linear', bounds_error=False, fill_value=0)
            if wlenhighres[0] < gcwlenmin:
                gcwlenmin = wlenhighres[0]
            if wlenhighres[-1] > gcwlenmax:
                gcwlenmax = wlenhighres[-1]
            carlastack = gcwlenintpol(wlenmultistack)
            multistack = multistack + carlastack

            #save the stack data:
            stackpath = 'stacking/figures/Stacking data/'
            run_name = 'testsaveall.fits'

            #check to see if file exists
            if os.path.isfile(stackpath+run_name): #read data if exists
                #open file to append new data for run
                hduold = fits.open(stackpath+run_name)
                oldcols = hduold[1].columns
                #prep new data
                runwlencol = fits.Column(name=gcname+'_Wlen', array=wlenmultistack, format='F')
                runnormspeccol = fits.Column(name=gcname+'_Stack_Flux', array=carlastack, format='F')
                newcols = fits.ColDefs([runwlencol, runnormspeccol])
                #combine and append to file
                runstackdata = fits.BinTableHDU.from_columns(oldcols+newcols)
                primary = fits.PrimaryHDU()
                hdul = fits.HDUList([primary, runstackdata])
                hduold.close()
                hdul.writeto(stackpath+run_name, overwrite = True)

            else: #only add data if not
                runwlencol = fits.Column(name=gcname+'_Wlen', array=wlenmultistack, format='F')
                runnormspeccol = fits.Column(name=gcname+'_Stack_Flux', array=carlastack, format='F')
                #combine and append to file
                runstackdata = fits.BinTableHDU.from_columns([runwlencol,runnormspeccol])
                primary = fits.PrimaryHDU()
                hdul = fits.HDUList([primary, runstackdata])
                hdul.writeto(stackpath+run_name, overwrite = True)


#cut down zero padding
gcstart = (np.abs(wlenmultistack - gcwlenmin)).argmin()
gcend = (np.abs(wlenmultistack - gcwlenmax)).argmin()
wlenmultistack = wlenmultistack[gcstart:gcend]
multistack = multistack[gcstart:gcend]

#plotting and data manipulation for output

meanstack = multistack/specstacktot
medianstack = multistack


#plot line in unstacked graph
plt.figure('multi-carla stack + uncombined')
plt.plot(np.array([(1215.67),(1215.67)]),np.array([150,-150]),'--')
plt.text(1215.67, 150, r' Ly$\alpha$ absorbtion for '+ str(carlatot) +' carla with '+ str(specstacktot) +' spectra')

plt.figure('stacked clusters')
plt.subplot(2,1,1)
plt.plot(wlenmultistack[:-300], meanstack[:-300])
plt.plot(np.array([(1215.67),(1215.67)]),np.array([np.min(meanstack),np.max(meanstack)]),'--')
plt.text(1215.67, np.min(meanstack), r' Ly$\alpha$ absorbtion for multi-carla stack with '+ str(specstacktot) +' spectra')
plt.xlabel(r' $\lambda$ ($\mathrm{\AA}$)')
plt.ylabel(r'$<F>$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')


plt.subplot(2,1,2)
plt.plot(wlenmultistack[:-300], medianstack[:-300])
plt.plot(np.array([(1215.67),(1215.67)]),np.array([np.min(meanstack),np.max(meanstack)]),'--')
plt.text(1215.67, np.min(meanstack), r' Ly$\alpha$ absorbtion for multi-carla stack with '+ str(specstacktot) +' spectra')
plt.xlabel(r' $\lambda$ ($\mathrm{\AA}$)')
plt.ylabel(r'$<F>$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')

plt.show()
