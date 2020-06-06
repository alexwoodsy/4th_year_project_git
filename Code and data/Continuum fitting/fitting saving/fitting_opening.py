
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
import os

#plt.style.use('mystyle')

specnames = next(os.walk('Fitted Spectra/'))[2]
spectot = len(specnames)
lya = 1215.67

for spec in specnames[1000:1100]:
    specdirectory = 'Fitted Spectra/' + spec
    data = fits.getdata(specdirectory, ext=1)
    metadata = fits.getdata(specdirectory, ext=2)
    redshift = metadata.field(0) #qso z

    wlen = data.field(0)
    flux = data.field(1)
    medcontfit = data.field(2)
    maxcontfit = data.field(3)

    wlenshift = wlen/(1+redshift)
    maxnormspec = (flux/maxcontfit)

    # #z =  (wlenshift/lya)*(1+redshift) - 1
    # z =  ((wlen/lya) -1)/(1+redshift) -1
    # a = 0.0018
    # b = 3.92
    # teff = a*(1+z)**b
    # faucher = np.exp(-teff)

    plate = spec[5:9]
    mjd = spec[10:15]
    fiberid = spec[16:20]

    path = 'E:/spectralyalpha/BOSSLyaDR9_spectra/BOSSLyaDR9_spectra/'+plate+'/'+'speclya-'+plate+'-'+mjd+'-'+fiberid+'.fits'

    hdul = fits.open(path)
    leedata = hdul[1].data
    contflag = hdul[1].header['CONTFLG']

    #leedata = fits.getdata(path,ext=1)#import fits image
    leeflux = leedata.field(0)
    leewlen = 10**leedata.field(1)
    leemodel = leedata.field(7)
    noisecorr = leedata.field(10)
    leecont = leedata.field(11)

    hdul.close()

    plt.plot(wlen,flux,label='ours')
    plt.plot(leewlen,leeflux,label='lee')
    plt.plot(wlen,maxcontfit,label='ours')
    plt.plot(leewlen,leecont,label='lee')

    plt.legend()



    plt.show()










#
