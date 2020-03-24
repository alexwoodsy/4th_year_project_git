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

#plt.style.use('mystyle') #path C:\Users\alexw\AppData\Local\Programs\Python\Python37\Lib\site-packages\matplotlib\mpl-data\stylelib

#imports the spectra from the spectra folder
specnames = next(os.walk('Spectra'))[2]
spectot = len(specnames)
#
# #set up sample to be looped over
# #read in data linking cluster to qso's
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
for carlaselect in range(0,1):
    specmatch = []
    for i in range(0,poslen): #CHNAGE
        if clusternames[i] == match[carlaselect] and specnames[i][-4:] == 'fits':
            specmatch.append(specnames[i])

    specmatch = ['spec-0343-51692-0145.fits','spec-0435-51882-0637.fits','spec-2947-54533-0417.fits','spec-3970-55591-0148.fits']
    normspeckstack = np.zeros(100000)
    wlenmin = []
    wlenmax = []

    gcredshift = 0#matchredshift[carlaselect]
    gclyalpha = 1215.67*(1+gcredshift)

    zlim = 2.2
    stonlim = 2

    wlentot = 0
    for spec in specmatch:
        spec = [spec]
        wlen, normspec, wlenlineind = fitmeth.contfitv7(spec, zlim , stonlim, showplot = True, showerror = True)
        wlenshift = wlen/(1+gcredshift)
        if len(wlenshift) > wlentot:
            wlentot = len(wlenshift)





    #     wlenmin.append(wlenshift[0])
    #     wlenmax.append(wlenshift[-1])
    # print('lim search done for ' + match[carlaselect])
    #
    # #zero pad spectra to the length of wlenhighres
    #
    # lowdiff, highdiff = wlenshift[0], wlenshift[-1]
    #
    #
    #
    # wlenhighres = np.linspace(np.max(wlenmin), np.min(wlenmax), 100000)
    #
    # for spec in specmatch:
    #     spec = [spec]
    #     wlen, normspec, lyalpha = fitmeth.contfitv7(spec, zlim ,stonlim, showplot = False, showerror = False)
    #     wlenshift = wlen/(1+gcredshift)
    #     wlenintpol = interpolate.interp1d(wlenshift, normspec, 'linear')
    #     normspechighres = wlenintpol(wlenhighres)
    #     plt.plot(wlenhighres,normspechighres)
    #     normspeckstack = normspeckstack + normspechighres
    # print('stacking done for ' + match[carlaselect])
    #
    #
    #
    # plt.figure('stacked spectra for '+ match[carlaselect])
    # plt.plot(wlenhighres, normspeckstack)
    # plt.show()
















#
