import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
from scipy.optimize import curve_fit as cf
import os, time
#import fitting
import sys

#plt.style.use('mystyle') #path C:\Users\alexw\AppData\Local\Programs\Python\Python37\Lib\site-packages\matplotlib\mpl-data\stylelib

def findval(array,val):
    array = np.asarray(array)
    ind = np.abs(array - val).argmin()
    return ind









#get lee anit-match from fits file made in topcat:
leefolderpath = 'E:/spectralyalpha/BOSSLyaDR9_spectra/BOSSLyaDR9_spectra/'

leeamdata = fits.getdata('Anti-Match/lee_anti-match.fits',ext=1)#import fits image
leeamlen = len(leeamdata)
leeam_z = []
leeam_sn = []



####################-----now stack them all FOR LEE------#####################
#runlen = leeamlen
runlen = 1000
#stacking initialisatin
fillval = np.nan
highreslen = 50000
cutinds = []
normspecstore = np.empty([runlen,highreslen])
contspecstore = np.empty([runlen,highreslen])
normspecstore[:] = fillval
contspecstore[:] = fillval

wlenhighres = np.linspace(500, 4500, highreslen)
wlenmin = 10000
wlenmax = 0
stackstatus = []
specnumber = 0


showerror = False

for ind in range(0,runlen):
    plate = str(leeamdata[ind][4]).zfill(4)
    mjd = str(leeamdata[ind][5]).zfill(5)
    fiberid = str(leeamdata[ind][6]).zfill(4)
    redshift = leeamdata[ind][7]
    stonall = leeamdata[ind][9]
    leeamfilename = 'speclya-'+plate+'-'+mjd+'-'+fiberid+'.fits'
    leeam_z.append(redshift)
    leeam_sn.append(stonall) #save for hist to show dist

    leeampath = leefolderpath+plate+'/'+leeamfilename
    #read in lee data:
    leedata = fits.getdata(leeampath,ext=1)#import fits image

    flux = leedata.field(0)
    wlen = 10**(leedata.field(1))
    contfit = leedata.field(11)

    lyalpha = (1+redshift)*1215.67
    lyalphaind = findval(wlen, lyalpha)
    gcredshift = 2.75 + np.random.random(1)-0.5
    gclyalpha = (1+gcredshift)*1215.67
    spec = leeamfilename

    zlims = np.array([gcredshift+0.05 , gcredshift + 2])
    stonlim = 5
    pw = 0


    #proceed with rest of the checks
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
        stackstatus.append('success')

        cut = tuple([contfit != 0 ])
        contfit = contfit[cut]
        flux = flux[cut]
        wlen = wlen[cut]

        wlenshift = wlen#/(1+gcredshift) #+np.random.random(1)-0.5
        normspec = flux/contfit

        # plt.plot(wlenshift,flux)
        # plt.plot(wlenshift,contfit)
        #plt.plot(wlenshift,normspec)
        #plt.show()

        wlenintpol = interpolate.interp1d(wlenshift, normspec, 'linear', bounds_error=False, fill_value=fillval)
        contintpol = interpolate.interp1d(wlenshift, contfit, 'linear', bounds_error=False, fill_value=fillval)
        if wlenshift[0] < wlenmin:
            wlenmin = wlenshift[0]
        if wlenshift[-1] > wlenmax:
            wlenmax = wlenshift[-1]

        normspecstore[specnumber, 0:] = wlenintpol(wlenhighres)
        contspecstore[specnumber, 0:] = contintpol(wlenhighres)

        specnumber = specnumber + 1


#output
print('stacking attempted for '+ str(len(stackstatus)) + ' spectra')

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

print(str(errortot) + ' did not meet conditions with:')
print(str(stacktot) +' stacked successfully with ' + str(errortot) + ' not used:')
print(str(zerrortot) + ' outside redshift range ' + str(zlims[0]) + ' < z < '+ str(zlims[1]))
print(str(foresterrortot) + ' lyalpaforest out of spectra range' )
print(str(gcerrortot) + ' galaxy luxter lyalpha out of spectra range' )
print(str(stonerrortot) + ' S/N below ' + str(stonlim))
print(' ')


#stacking data only if there are spec to stack

#cut extra zeropadding
start = (np.abs(wlenhighres - wlenmin)).argmin() +10
end = (np.abs(wlenhighres - wlenmax)).argmin() -10
wlenhighres = wlenhighres[start:end]
#cut downcols / remove empty rows
#normspecstore = np.delete(normspecstore, cutinds, axis = 0)
normspecstore = normspecstore[0:, start:end]
meanspec = np.nanmean(normspecstore, axis=0)
medspec = np.nanmedian(normspecstore, axis=0)

contspecstore = contspecstore[0:, start:end]

meancont = np.nanmean(contspecstore, axis=0)
medcont = np.nanmedian(contspecstore, axis=0)

plt.plot(wlenhighres,meanspec)
plt.show()









#
