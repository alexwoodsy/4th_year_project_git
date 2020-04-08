import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate
import os

#import metadata tableand data
metatabledata = fits.getdata('metatable.fits', ext=1)
fitlen = len(metatabledata)

metaspecfilename = []
redshift = np.zeros(fitlen)
stonall = np.zeros(fitlen)
stonforest = np.zeros(fitlen)
lyalpha = np.zeros(fitlen)
lyalphaind = np.zeros(fitlen)
metagcname = []
gcredshift = np.zeros(fitlen)
metagc_qso_sep = np.zeros(fitlen)
gclyalpha = np.zeros(fitlen)
gclyalphaind = np.zeros(fitlen)
stackmsg = []


for i in range(0,fitlen):
    #extract metadata
    metaspecfilename.append(metatabledata[i][0])
    redshift[i]= metatabledata[i][1]
    stonall[i] = metatabledata[i][2]
    stonforest[i] = metatabledata[i][3]
    lyalpha[i] = metatabledata[i][4]
    lyalphaind[i] = metatabledata[i][5]
    metagcname.append(metatabledata[i][6])
    gcredshift[i] = metatabledata[i][7]
    metagc_qso_sep[i] = metatabledata[i][8]
    gclyalpha[i] = metatabledata[i][9]
    gclyalphaind[i] = metatabledata[i][10]
    stackmsg.append(metatabledata[i][11])

#get carla sample
carlamatchnew = list(set(metagcname))
print(carlamatchnew)
print(len(carlamatchnew))

###########old method ##############

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
    if c > 0:
        carlamatch.append(carlanames[i])
print(carlamatch)
print(len(carlamatch))

matches = []
for specnew in carlamatchnew:
    count = 0
    for specold in carlamatch:
        if specold == specnew:
            count = 'match'
    if count != 'match':
        matches.append(specnew)

print(matches)
print(len(matches))
