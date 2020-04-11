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
overdensity = np.zeros(carlalen)
carlamatch =  []
for i in range(0,carlalen):
    c = 0
    carlanames = str(carladata[i][0])
    overdensity[i] = carladata[i][6]
    for k in range(0,fitlen):
        if carlanames == metagcname[k]:# and overdensity[10]== overdensity[i]: #add conditional for overdendity etc.. here
            c = c + 1
    if c > 0:
        carlamatch.append(carlanames)


specstacktot = 0 #total number of spectra stacked in all carla

carlarange = np.arange(0,len(carlamatch))
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
            zlims = np.array([gcredshift + 0.1, gcredshift + 5])
            stonlim = 1
            pw = 0
            #check fit can be added to stack:
            if gc_qso_sep > 240: #ignore qso if it has no forest (lyalpha - wlen[0] <= 0)
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
                wlenintpol = interpolate.interp1d(wlenshift, normspec, 'linear', bounds_error=False, fill_value=0)
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
            meangcwlenintpol = interpolate.interp1d(wlenhighres, meanspec, 'linear', bounds_error=False, fill_value=0)
            medgcwlenintpol = interpolate.interp1d(wlenhighres, medspec, 'linear', bounds_error=False, fill_value=0)
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


#plt.figure('('+str(carlanumber)+'/'+str(carlafail)+ ') stacked clusters')
fig1, ax = plt.subplots(2,1,num='multi-carla stack - uncombined')
for c in range(0,carlasucces):
    ax[0].plot(wlenmultistack[:-300], meanmultistore[c, :-300])
    ax[0].set_xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
    ax[0].set_ylabel(r'$<F>$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')

    ax[1].plot(wlenmultistack[:-300], medmultistore[c, :-300])
    ax[1].set_xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
    ax[1].set_ylabel(r'MEDIAN $F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')



fig2, ax = plt.subplots(2,1)
fig2.canvas.set_window_title('('+str(carlasucces)+'/'+str(carlanumber)+ ') stacked clusters')
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
