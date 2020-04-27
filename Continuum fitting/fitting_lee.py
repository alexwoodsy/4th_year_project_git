import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate
import os


folderpath = 'E:/spectralyalpha/BOSSLyaDR9_spectra/BOSSLyaDR9_spectra/'

matchdata = fits.getdata('matchsmaple.fits',ext=1)#import fits image
matchlen = len(matchdata)
speclyapath = []
specpath = []

for i in range(0,matchlen):
    plate = str(matchdata[i][4]).zfill(4)
    mjd = str(matchdata[i][5]).zfill(5)
    fiberid = str(matchdata[i][6]).zfill(4)
    path = folderpath+plate+'/'+'speclya-'+plate+'-'+mjd+'-'+fiberid+'.fits'
    speclyapath.append(path)

    #get corresponding name in postable (dr14)
    newplate = str(matchdata[i][20]).zfill(4)
    newmjd = str(matchdata[i][21]).zfill(5)
    newfiberid = str(matchdata[i][22]).zfill(4)
    newpath = 'spec-'+newplate+'-'+newmjd+'-'+newfiberid+'.fits'
    specpath.append(newpath)

specmatchlya =[]

for j in range(0,len(specmatch)):
    specfilename = specmatch[j]
    #finds matching spec
    for i in range(0,matchlen):
        check = specpath[i]
        if check == specfilename:
            out = speclyapath[i]
    specmatchlya.append(out)

print(len(specmatchlya))    


#
# #show matching spec
# for i in range(len(specmatchlya)):
#     if specmatchlya[i] != 'null':
#         data = fits.getdata(specmatchlya[i],ext=1)#import fits image
#         speclen = len(data)
#         flux = np.zeros(speclen)
#         wlen = np.zeros(speclen)
#         model = np.zeros(speclen)
#         ivar = np.zeros(speclen)
#         cont = np.zeros(speclen)
#
#         for i in range(0,speclen):
#             flux[i] = data[i][0]
#             ivar[i] = data[i][2]
#             if ivar[i] == 0:
#                 ivar[i] = 0.00001
#             wlen[i] = 10**(data[i][1])
#             model[i] = data[i][7]
#             cont[i] = data[i][11]
#
#         contrem =flux-cont
#
#         plt.figure('flux and fitting')
#         plt.plot(wlen,flux, label='flux')
#         plt.plot(wlen,cont, label='cont fit')
#         plt.legend()
#
#         plt.figure('continuum subtraction')
#         plt.plot(wlen,contrem, label='continuum removed')
#         plt.legend()
#         plt.show()
#
#
#









 #
