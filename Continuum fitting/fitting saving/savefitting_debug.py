import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
import os
#import fitting
import sys
#sys.path.append('C:/Users/jason/GIT/4th_year_project_git/Continuum Fitting')
sys.path.append('C:/Users/alexw/Documents/GitHub/4th_year_project_git/Continuum fitting')


#imports the spectra from the spectra folder
specnames = next(os.walk('Fitted Spectra'))[2]
spectot = len(specnames)

oldspecnames = next(os.walk('Spectra'))[2]
oldspectot = len(specnames)


nullerror = [] ###135 null error
count = 0
for spec in specnames:
    if len(spec) != 35: #standard length == 35 hence anything more == spec in 2 clusters
        print(spec[21:-16])
        
    fitdata = fits.getdata('Fitted Spectra/'+spec, ext=2)
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

#     if gcname == 'null':
#         count = count + 1
#         #print('found ' + str(count))
#         nullerror.append(spec)
#
# print(count)
# print(nullerror)

#delete old null error files












######carla stuff
#set up sample to be looped over
#read in data linking cluster to qso's
# posdata = fits.getdata('PositionsTable.fits',ext=1)
# poslen = len(posdata)
#
# clusternames = np.zeros(poslen).astype(str)
# specfilename = np.zeros(poslen).astype(str)
# specfilenamefitted = np.zeros(poslen).astype(str)
# for j in range(0,poslen):
#     #get cluster info
#     clusternames[j] = str(posdata[j][105])
#     #get assocaited id to find spec in folder
#     plate = str(posdata[j][4]).zfill(4)
#     mjd = str(posdata[j][5]).zfill(5)
#     fiberid = str(posdata[j][6]).zfill(4)
#     specfilename[j] = 'spec-'+plate+'-'+mjd+'-'+fiberid+'.fits'
#     specfilenamefitted[j] = 'spec-'+plate+'-'+mjd+'-'+fiberid
#
# #read in carla
# carladata = fits.getdata('CARLA/CARLA_table_Aug_2013.fits',ext=1)
# carlalen = len(carladata)
# carlanames = np.zeros(carlalen).astype(str)
#
# for i in range(0,carlalen):
#     carlanames[i] = str(carladata[i][0])
#
#
# #select carla agn that were selected in filtering (in pos table) and associated spec index in table
# carlamatch =  []
#
# for i in range(0,carlalen):
#     c = 0
#     for k in range(0,poslen):
#         if carlanames[i] == clusternames[k]:
#             c = c + 1
#     if c > 1:
#         carlamatch.append(carlanames[i])
#
# carlaselect = 0
#
#
# specmatchfitted =[]
# for i in range(0,poslen): #CHNAGE
#     if clusternames[i] == carlamatch[carlaselect]:
#         output = specfilenamefitted[i]+'-prefitted.fits'
#         specmatchfitted.append(output)
################

#
# count = 0
# notfound = []
#
# for specold in oldspecnames:
#     count = 0
#     for spec in specnames:
#         if specold[:-5] == spec[:-15]:
#             count = count + 1
#     if count == 0:
#         notfound.append(specold)
#
# ###185 not in spec folder
# print(notfound)



#
