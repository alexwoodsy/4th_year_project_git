#stacking v2
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
import os
#import fitting
import sys
sys.path.append('C:/Users/jason/GIT/4th_year_project_git/Continuum Fitting')

# sys.path.append('C:/Users/alexw/Documents/GitHub/4th_year_project_git/Continuum fitting')
#path for other pc sys.path.append('C:/Users/alexw/OneDrive/Documents/University work/4th year work/Main project/4th_year_project_git/Continuum fitting')
import fitting_jason as fitmeth

plt.style.use('mystyle')

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
clusterstack = []
clusterhighres = np.linspace(0, 6000, 1000000)
clustermin = 10000
clustermax = 0
#get spec names for a carla target - do stacking
for carlaselect in range(0,1):
    specmatch = []
    for i in range(0,poslen): #CHNAGE
        if clusternames[i] == match[carlaselect] and specnames[i][-4:] == 'fits':
            specmatch.append(specnames[i])

    gcredshift = matchredshift[carlaselect]

    zlim = gcredshift+0.15
    stonlim = 1

    specmatch = specmatch#[0:20]
    #specmatch = ['spec-0343-51692-0145.fits','spec-0435-51882-0637.fits','spec-2947-54533-0417.fits','spec-3970-55591-0148.fits']

    normspeckstack = []
    wlenhighres = np.linspace(0, 6000, 1000000)
    wlenmin = 10000
    wlenmax = 0
    for spec in specmatch:
        spec = [spec]
        wlen, normspec, wlenlineind, redshift,stackstatus = fitmeth.contfitv7(spec, zlim , stonlim, gcredshift, showplot = False, showerror = False)
        wlenshift = wlen/(1+gcredshift)
        wlenintpol = interpolate.interp1d(wlenshift, normspec, 'linear', bounds_error=False, fill_value=0)
        if wlenshift[0] < wlenmin:
            wlenmin = wlenshift[0]
        if wlenshift[-1] > wlenmax:
            wlenmax = wlenshift[-1]

        normspechighres = wlenintpol(wlenhighres)
        # print(normspeckstack)
        normspeckstack = np.append(normspeckstack,normspechighres,axis=0)
        # normspecmed.append(normspechighres)
        # print(np.max(normspeckstack))
        # plt.figure('Unstacked Spectra')
        plt.plot(wlenshift, normspec)
        plt.xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
        plt.ylabel(r'$F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')
    print('stacking done for ' + match[carlaselect])
    clusterstack = np.append(clusterstack,normspeckstack,axis=0)
    # clustermed.append(normspecmed)
    if wlenmin < clustermin:
        clustermin = wlenmin
    if wlenmax > clustermax:
        clustermax = wlenmax
    # normspeckavg = normspeckstack/len(specmatch)
    # cut extra zeropadding
# print(clusterstack)
clustermed = np.median(normspeckstack,axis=0)
print(clustermed)

clustermax = 1600
# clusterstack = np.median(clustermed)
start = (np.abs(clusterhighres - clustermin)).argmin()
end = (np.abs(clusterhighres - clustermax)).argmin()
clusterhighres = clusterhighres[start:end]
clustermed = clustermed[start:end]

#downsample stacked specind
# dsrange = np.linspace(clusterhighres[0], clusterhighres[-1],5000)
dsrange = np.linspace(clusterhighres[0],clusterhighres[-1],5000)
dsnormspewcstack = signal.resample(clustermed, 5000)

plt.figure('Stacked Spectra')
plt.plot(dsrange, dsnormspewcstack)
plt.plot([1216.57,1216.57],[np.min(dsnormspewcstack)*1.1,np.max(dsnormspewcstack)*1.1],'r--')
plt.xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
plt.ylabel(r'$F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$')
plt.text(1215,np.max(dsnormspewcstack)*1.11,r'ly$\alpha$',fontsize = 10)
plt.show()
