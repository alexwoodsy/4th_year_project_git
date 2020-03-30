#stacking v2
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
import os
#import fitting
import sys
sys.path.append('C:/Users/alexw/Documents/GitHub/4th_year_project_git/Continuum fitting')
#path for other pc sys.path.append('C:/Users/alexw/OneDrive/Documents/University work/4th year work/Main project/4th_year_project_git/Continuum fitting')
import fitting_v7 as fitmeth

plt.style.use('mystyle') #path C:\Users\alexw\AppData\Local\Programs\Python\Python37\Lib\site-packages\matplotlib\mpl-data\stylelib

#imports the spectra from the spectra folder
specnames = next(os.walk('Spectra'))[2]
spectot = len(specnames)

#set up sample to be looped over
#read in data linking cluster to qso's
posdata = fits.getdata('PositionsTable.fits',ext=1)
poslen = len(posdata)
clusterredshift = np.zeros(poslen)
clusternames = np.zeros(poslen).astype(str)
specfilename = np.zeros(poslen).astype(str)

for j in range(0,poslen):
    #get cluster info
    clusterredshift[j] = posdata[j][109]
    clusternames[j] = str(posdata[j][105])
    #get assocaited id to find spec in folder
    plate = str(posdata[j][4]).zfill(4)
    mjd = str(posdata[j][5]).zfill(5)
    fiberid = str(posdata[j][6]).zfill(4)
    specfilename[j] = 'spec-'+plate+'-'+mjd+'-'+fiberid+'.fits'

#read in carla
carladata = fits.getdata('CARLA/CARLA_table_Aug_2013.fits',ext=1)
carlalen = len(carladata)
carlanames = np.zeros(carlalen).astype(str)
carlaredshift = np.zeros(carlalen)
for i in range(0,carlalen):
    carlanames[i] = str(carladata[i][0])
    carlaredshift[i] = carladata[i][4]


#select carla agn that were selected in filtering (in pos table) and associated spec index in table
match =  []
matchredshift =[]

for i in range(0,carlalen):
    c = 0
    for k in range(0,poslen):
        if carlanames[i] == clusternames[k]:
            c = c + 1
    if c > 1:
        match.append(carlanames[i])
        matchredshift.append(carlaredshift[i])


#####----STACKING---#####

#get spec names for a carla target - do stacking
for carlaselect in range(0,1): #all change to - matchlen
    print('processing '+ match[carlaselect])

    specmatch = []
    for i in range(0,poslen): #CHNAGE
        if clusternames[i] == match[carlaselect] and specnames[i][-4:] == 'fits':
            specmatch.append(specnames[i])


    #gc absorber location
    gcredshift = matchredshift[carlaselect]

    zlim = 2.15
    stonlim = 1

    #specmatch = specmatch[0:5]
    specmatch = ['spec-0343-51692-0145.fits','spec-0435-51882-0637.fits','spec-2947-54533-0417.fits','spec-3970-55591-0148.fits']

    normspeckstack = np.zeros(1000000)
    wlenhighres = np.linspace(0, 6000, 1000000)
    wlenmin = 10000
    wlenmax = 0
    stackstatus = []
    for spec in specmatch:
        spec = [spec]
        wlen, normspec, wlenlineind, redshift, stackcode = fitmeth.contfitv7(spec, zlim , stonlim, gcredshift,showplot = True, showerror = True)
        stackstatus.append(stackcode)
        wlenshift = wlen/(1+gcredshift)
        wlenintpol = interpolate.interp1d(wlenshift, normspec, 'linear', bounds_error=False, fill_value=0)
        if wlenshift[0] < wlenmin:
            wlenmin = wlenshift[0]
        if wlenshift[-1] > wlenmax:
            wlenmax = wlenshift[-1]
        normspechighres = wlenintpol(wlenhighres)
        normspeckstack = normspeckstack + normspechighres
        if np.mean(normspec) != 0:
            plt.figure(match[carlaselect] +' unstacked ')
            plt.plot(wlenshift, normspec)
            plt.xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
            plt.ylabel(r'$F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')

    figuretitle = match[carlaselect] + '_unstacked'
    savepath = 'stacking/figures/'+figuretitle+'.svg'
    plt.savefig(savepath)

    #output
    print(' ')
    print('------------------------------------------------------------------------------------------------------------')
    print(' ')
    print('stacking attempted for '+ str(len(stackstatus)) + ' spectra for CARLA: ' + match[carlaselect])

    stacktot = zerrortot = foresterrortot = stonerrortot = 0
    for x in stackstatus:
        if x == 'success':
            stacktot = stacktot + 1
        if x == 'zerror':
            zerrortot = zerrortot + 1
        if x == 'foresterror':
            foresterrortot = foresterrortot + 1
        if x =='stonerror':
            stonerrortot = stonerrortot + 1
    errortot = zerrortot + foresterrortot + stonerrortot

    print(str(stacktot) +' stacked successfully with ' + str(errortot) + ' not used with:')
    print(str(zerrortot) + ' outside redshift range z > ' + str(zlim))
    print(str(foresterrortot) + ' lyalpaforest out of spectra range' )
    print(str(stonerrortot) + ' S/N below ' + str(stonlim))


    #stacking data:

    #cut extra zeropadding
    start = (np.abs(wlenhighres - wlenmin)).argmin()
    end = (np.abs(wlenhighres - wlenmax)).argmin()
    wlenhighres = wlenhighres[start:end]
    normspeckstack = normspeckstack[start:end]

    #downsample stacked specind
    dsrange = np.linspace(wlenhighres[0], wlenhighres[-1],5000)
    dsnormspewcstack = signal.resample(normspeckstack, 5000)

    figuretitle = match[carlaselect] + '_stacked'
    plt.figure(figuretitle)
    plt.plot(dsrange, dsnormspewcstack)
    plt.xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
    plt.ylabel(r'$F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')
    string = match[carlaselect] + ' z = '+ str(gcredshift) + r' Ly$\alpha$ absorbtion'
    plt.text((1215.67), np.max(dsnormspewcstack), string)
    plt.plot(np.array([(1215.67),(1215.67)]),np.array([np.min(dsnormspewcstack),np.max(dsnormspewcstack)]),'--')
    savepath = 'stacking/figures/'+figuretitle+'.svg'
    #plt.savefig(savepath)
    #plt.close('all')
    plt.show()
    print('saving '+ figuretitle)










    #
