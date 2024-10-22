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


######carla stuff
#set up sample to be looped over
#read in data linking cluster to qso's
posdata = fits.getdata('PositionsTable.fits',ext=1)
poslen = len(posdata)

clusternames = np.zeros(poslen).astype(str)
specfilename = np.zeros(poslen).astype(str)
specfilenamefitted = np.zeros(poslen).astype(str)
for j in range(0,poslen):
    #get cluster info
    clusternames[j] = str(posdata[j][105])
    #get assocaited id to find spec in folder
    plate = str(posdata[j][4]).zfill(4)
    mjd = str(posdata[j][5]).zfill(5)
    fiberid = str(posdata[j][6]).zfill(4)
    specfilename[j] = 'spec-'+plate+'-'+mjd+'-'+fiberid+'.fits'
    specfilenamefitted[j] = 'spec-'+plate+'-'+mjd+'-'+fiberid

#read in carla
carladata = fits.getdata('CARLA/CARLA_table_Aug_2013.fits',ext=1)
carlalen = len(carladata)
carlanames = np.zeros(carlalen).astype(str)

for i in range(0,carlalen):
    carlanames[i] = str(carladata[i][0])


#select carla agn that were selected in filtering (in pos table) and associated spec index in table
carlamatch =  []

for i in range(0,carlalen):
    c = 0
    for k in range(0,poslen):
        if carlanames[i] == clusternames[k]:
            c = c + 1
    if c > 1:
        carlamatch.append(carlanames[i])


#####----STACKING---#####

#get spec names for a carla target - do stacking

wlenmultistack = np.linspace(0, 6000, 100000)
multistack = np.zeros(100000)
gcwlenmin = 10000
gcwlenmax = 0


specstacktot = 0

for carlaselect in range(0,1): #all change to - matchlen
    print('processing '+ carlamatch[carlaselect])
    #get carla spectra to loop over
    specmatchfitted =[]
    for i in range(0,poslen): #CHNAGE
        if clusternames[i] == carlamatch[carlaselect]:
            output = specfilenamefitted[i]+'-prefitted.fits'
            specmatchfitted.append(output)

    normspeckstack = np.zeros(50000)
    wlenhighres = np.linspace(0, 6000, 50000)
    wlenmin = 10000
    wlenmax = 0
    stackstatus = []

    #debugging parameters:
    showerror = True
    showplot = False
    for spec in specmatchfitted:
        specdirectory = 'Fitted Spectra/'+spec
        data = fits.getdata(specdirectory,ext=1)
        speclen = len(data)
        fitdata = fits.getdata(specdirectory,ext=2)#read in fits meta data
        #predefine data variables
        normspec = np.zeros(speclen)
        wlen = np.zeros(speclen)
        #extract fitting data
        for i in range(0,speclen):
            normspec[i] = data[i][1]
            wlen[i] = (data[i][0])
        #extract metadata
        redshift= fitdata[0][0]
        stonall = fitdata[0][1]
        stonforest = fitdata[0][2]
        lyalpha = fitdata[0][3]
        lyalphaind = fitdata[0][4]
        gcname = fitdata[0][5]
        gcredshift = fitdata[0][6]
        gc_qso_sep = fitdata[0][7]
        gclyalpha = fitdata[0][8]
        gclyalphaind = fitdata[0][9]
        stackmsg = fitdata[0][10]
        print(spec +' redshift = '+ str(redshift) +' gc = '+ gcname)

        #conditions for stacking
        zlims = np.array([gcredshift + 0.5, gcredshift + 2])
        stonlim = 1
        pw = 0

        #check fit can be added to stack:
        if lyalpha - wlen[0] <= 0: #ignore qso if it has no forest (lyalpha - wlen[0] <= 0)
            stackcode = 'foresterror'
            if showerror == True:
                print('forest warning!: '+spec+' z = '+ str(redshift) +' forest before start of spectrum (S/N = '+ str(stonall) +')')
        elif redshift < zlims[0] and redshift > zlims[1]: #removes qso in carla (those within z = 0.05) and too high
            stackcode = 'zerror'
            if showerror == True:
                print('zlim warning!: '+spec+' z = '+ str(redshift) +' below not in limits (S/N = '+ str(stonall) +')')
        elif stonall < stonlim:
            stackcode = 'stonerror'
            if showerror == True:
                print('S/N warning!: '+spec+' S/N = '+ str(stonall) +' too low (z = '+ str(redshift) +')')
        elif gclyalpha - wlen[0] <= pw:
            stackcode = 'gcerror'
            if showerror == True:
                print('gc warning!: '+spec+' has z too low with respect to gc withlyalpha '+ str(gclyalpha - wlen[0]) + ' Angstroms before start of spectra')
        else:
            #proceed with continuum fitting
            stackcode = 'success'
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
            figuretitleunstacked = carlamatch[carlaselect] + '_unstacked'
            plt.figure(figuretitleunstacked)
            plt.plot(wlenshift, normspec)
            plt.xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
            plt.ylabel(r'$F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')

    stackstatus.append(stackcode)
    #savepath = 'stacking/figures/'+figuretitleunstacked+'.svg'
    #plt.savefig(savepath)
    #plt.close(figuretitleunstacked)

    #output
    print(' ')
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

    print(str(stacktot) +' stacked successfully with ' + str(errortot) + ' not used:')
    print(str(zerrortot) + ' outside redshift range ' + str(zlims[0]) + ' < z < '+ str(zlims[1]))
    print(str(foresterrortot) + ' lyalpaforest out of spectra range' )
    print(str(gcerrortot) + ' galaxy luxter lyalpha out of spectra range' )
    print(str(stonerrortot) + ' S/N below ' + str(stonlim))
    specstacktot = specstacktot + stacktot

    #stacking data:
    #cut extra zeropadding
    start = (np.abs(wlenhighres - wlenmin)).argmin()
    end = (np.abs(wlenhighres - wlenmax)).argmin()
    wlenhighres = wlenhighres[start:end]
    normspeckstack = normspeckstack[start:end]


    # figuretitle = match[carlaselect] + '_stacked'
    # plt.figure(figuretitle)
    # plt.plot(wlenhighres, normspeckstack)
    # plt.xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
    # plt.ylabel(r'$F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')
    # string = match[carlaselect] + ' z = '+ str(gcredshift) + r' Ly$\alpha$ absorbtion'
    # plt.text((1215.67), np.max(normspeckstack), string)
    # plt.plot(np.array([(1215.67),(1215.67)]),np.array([np.min(normspeckstack),np.max(normspeckstack)]),'--')

    #savepath = 'stacking/figures/'+figuretitle+'.svg'
    #plt.savefig(savepath)
    #print('saving '+ figuretitle)
    #plt.close(figuretitle)

    #dsrange = np.linspace(wlenhighres[0], wlenhighres[-1],5000)
    #dsnormspewcstack = signal.resample(normspeckstack, 5000)

    plt.figure('multi-carla stack + uncombined')
    plt.subplot(2, 1, 1)
    plt.plot(wlenhighres, normspeckstack, label = carlamatch[carlaselect] + ' stack of '+str(stacktot) + ' QSO')
    plt.xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
    plt.ylabel(r'$F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')
    plt.legend()

#stack carla
gcwlenintpol = interpolate.interp1d(wlenhighres, normspeckstack, 'linear', bounds_error=False, fill_value=0)
if wlenhighres[0] < gcwlenmin:
    gcwlenmin = wlenhighres[0]
if wlenhighres[-1] > gcwlenmax:
    gcwlenmax = wlenhighres[-1]
carlastack = gcwlenintpol(wlenmultistack)
multistack = multistack + carlastack
#cut down zero padding
gcstart = (np.abs(wlenmultistack - gcwlenmin)).argmin()
gcend = (np.abs(wlenmultistack - gcwlenmax)).argmin()
wlenmultistack = wlenmultistack[gcstart:gcend]
multistack = multistack[gcstart:gcend]

pltrange = np.arange(0,int(len(multistack)*1))
meanstack = (multistack[pltrange]/1)
#normstack = meanstack/np.max(abs(meanstack))

plt.subplot(2, 1, 2)
plt.plot(wlenmultistack[pltrange], meanstack[pltrange])
plt.plot(np.array([(1215.67),(1215.67)]),np.array([np.min(meanstack),np.max(meanstack)]),'--')
string = r' Ly$\alpha$ absorbtion for multi-carla stack with '+ str(specstacktot) +' spectra'
plt.text(1215.67, 5, string)
plt.xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
plt.ylabel(r'$<F>$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')


plt.show()
