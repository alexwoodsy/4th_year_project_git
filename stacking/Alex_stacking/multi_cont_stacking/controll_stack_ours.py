import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
from scipy.optimize import curve_fit as cf
import os, random
#import fitting
import sys

#plt.style.use('mystyle') #path C:\Users\alexw\AppData\Local\Programs\Python\Python37\Lib\site-packages\matplotlib\mpl-data\stylelib

def findval(array,val):
    array = np.asarray(array)
    ind = np.abs(array - val).argmin()
    return ind



#################debugging selection parameters parameters:########################################
##########-----input stack -------###################################################
runsavename = 'ours_0to1000bin'
########---show stacking outpit-------#####################################################
showerror = False #see continuum fits + inividiual spec stack info
#################------save controll stack data:-----------########################################
saveoutput = True
#######-amount of spectra to attempt to stack based on input distro-###############################
binmult = 10 #multiplier of spec in each input stack bin to take from the anitmatch same
###################################################################################################


#read in input stackmetadata:
inname = 'stacking/figures/Stacking data/' + runsavename+ '.fits'
carladata = fits.getdata(inname,ext=2)
carlanames = carladata.field(0)
carlanumber = len(carladata)
qsodata = fits.getdata(inname,ext=3)
qsonumber = len(qsodata)
qso_z = qsodata.field(1)
qso_ston = qsodata.field(2)
qso_gcz = qsodata.field(3)
qso_gcname = qsodata.field(4)


#read in antimatchfiles metadata
folderpath = 'E:/AM-Prefitted Spectra/'#BIG PC ROUTE
#folderpath = 'C:/Users/alexw/Spectra_files/AM-Prefitted Spectra/AM-Prefitted Spectra/' #surface route

amdata = fits.getdata('Anti-Match/AM-PLUS_SN.fits',ext=1)#import fits image
amlen = len(amdata)#
am_names = amdata.field(0)
am_plate = amdata.field(3)
am_mjd = amdata.field(4)
am_fiberid = amdata.field(5)
am_z = amdata.field(6)
am_stonall = amdata.field(7)
am_stonforest = amdata.field(8)

#distribution information

qso_zdist = qso_z
zmin = np.floor(10*np.min(qso_zdist))/10
zmax = np.ceil(10*np.max(qso_zdist))/10

binwidth = 0.05
bins = np.arange(zmin,zmax,binwidth)
zbinsinds = np.digitize(qso_zdist, bins) - 1
stack_zbincounts = np.bincount(zbinsinds, minlength=len(bins))

plt.figure('norm zdists')
plt.step(bins+binwidth,stack_zbincounts/np.max(stack_zbincounts),label='stack dist')

#bin gc z info:
gcbins = np.arange(zmin,zmax,binwidth)
gczbinsinds = np.digitize(qso_gcz, gcbins)
stack_gczbincounts = np.bincount(gczbinsinds, minlength=len(gcbins))
plt.step(gcbins+binwidth,stack_gczbincounts/np.max(stack_gczbincounts),label='stack gc dist')

#bin the am spec in the same way:
am_zorderind = np.argsort(am_z).astype(int)
amz_start = findval(am_z,zmin) # get cut val
amz_end = findval(am_z,zmax) # get cut val
am_zorderind = am_zorderind[findval(am_zorderind,amz_start):findval(am_zorderind,amz_end)] #find cut val in orderd inds then trim the ordered inds and apply:

am_zcutdown = am_z[am_zorderind] #apply the cut
am_zbinsinds = np.digitize(am_zcutdown, bins) - 1 #bin each value in am_zcutdown belongs too
am_zbincounts = np.bincount(am_zbinsinds, minlength=len(bins))
print(bins)
print(am_zbincounts)

plt.step(bins+binwidth,am_zbincounts/np.max(am_zbincounts),label='amz dist')
plt.legend()
#plt.show()


print('cut am mean = '+str(np.mean(am_zcutdown))+' and std = '+ str(np.std(am_zcutdown)))
print('stack mean = '+str(np.mean(qso_zdist))+' and std = '+str(np.std(qso_zdist)))
print('stack gc mean = '+str(np.mean(qso_gcz))+' and std = '+str(np.std(qso_gcz)))
#print('am mean = '+str(np.mean(am_z))+' and std = '+ str(np.std(am_z)))
print('|')




####################-----now stack them all FOR OURS------#####################
fillval = np.nan
carlarange = np.arange(0,carlanumber)
#variables for stacking carla together
carlahighreslen = 100000
wlenmultistack = np.linspace(500, 4000, carlahighreslen)
carlacutinds = []
meanmultistore = np.empty([len(carlarange),carlahighreslen])
medmultistore = np.empty([len(carlarange),carlahighreslen])
meanmultistore[:] = fillval
medmultistore[:] = fillval

gcwlenmin = carlahighreslen
gcwlenmax = 0

stackedinds = [] #append inds of stacked am to ensure no double stacking
qsostacked = []
qsostacked_z = []
qsostacked_sn = []


for carla in carlarange:
    print('processing antimatch qso for gc in stack '+ carlanames[carla] + ' ('+str(carla)+'/'+str(carlanumber)+')')
    ############-----determine qso from antimtach to stack------###################
    cinds = np.array([]).astype(int)
    for ind in range(0, qsonumber):
        if qso_gcname[ind] == carlanames[carla]:
            cinds = np.append(cinds,ind)
    qso_zdist = qso_z[cinds]
    gcredshift = qso_gcz[cinds[0]]
    #print(gcredshift)

    #qso_gczdist = qso_gcz[cinds] all the same as same carla


    zmin = np.floor(10*np.min(qso_zdist))/10
    zmax = np.ceil(10*np.max(qso_zdist))/10

    binwidth = 0.05
    scbins = np.arange(zmin,zmax,binwidth)
    #print(scbins)
    sczbinsinds = np.digitize(qso_zdist, scbins) - 1
    sc_zbincounts = np.bincount(sczbinsinds, minlength=len(scbins))


    ########-----find matched in am by z and map matching gc deredshift value to am-spec------#######

    #inds of qso in amdata to stack defined by input stack metadata
    specstackinds = []
    for binind in range(0,len(bins)):
        #print('sc '+ str(round(bins[binind],2)))
        for scbinind in range(0,len(scbins)):
            #print(round(scbins[scbinind],2))
            if round(bins[binind],2) == round(scbins[scbinind],2):
                bin_scnumber = sc_zbincounts[scbinind]
                #print(scbins[scbinind])
                bin_aminds = list(am_zorderind[am_zbinsinds == binind]) #am spec z index's in the bins of the carla we are looking at
                #print(len(bin_aminds))
                if bin_scnumber !=0 and len(bin_aminds) !=0:
                    #take an amount in each bin given by input stack distribution and multiplier
                    choicenumber = bin_scnumber*binmult
                    bin_selection = random.sample(bin_aminds,k=choicenumber)
                    #print(bin_selection)
                    #append to list to stack
                    specstackinds = specstackinds + bin_selection

    #append bin selection to index store so they arent used in next carla
    #stackedinds = stackedinds + specstackinds
    #print(len(specstackinds))


    ######-----stacking code----#######
    runlen = len(specstackinds)
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


    for ind in specstackinds:
        plate = str(am_plate[ind]).zfill(4)
        mjd = str(am_mjd[ind]).zfill(5)
        fiberid = str(am_fiberid[ind]).zfill(4)
        redshift = am_z[ind]
        stonall = am_stonall[ind]
        stonforest = am_stonall[ind]
        spec = 'AM-spec-'+plate+'-'+mjd+'-'+fiberid+'-Preffited.fits'

        ampath = folderpath+spec
        #read in lee data:
        specdata = fits.getdata(ampath,ext=1)#import fits image

        wlen = specdata.field(0)
        flux = specdata.field(1)
        contfit = specdata.field(2)

        lyalpha = (1+redshift)*1215.67
        lyalphaind = findval(wlen, lyalpha)
        #gcredshift = 2.35 + np.random.random(1)-0.5
        gclyalpha = (1+gcredshift)*1215.67

        zlims = np.array([gcredshift+0.05 , gcredshift + 2])
        stonlim = 4
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
        elif stonforest < stonlim:
            stackstatus.append('stonerror')
            if showerror == True:
                print('S/N warning!: '+spec+' S/N = '+ str(stonforest) +' too low (z = '+ str(redshift) +')')
        elif gclyalpha - wlen[0] <= pw:
            stackstatus.append('gcerror')
            if showerror == True:
                print('gc warning!: '+spec+' has z too low with respect to gc withlyalpha '+ str(gclyalpha - wlen[0]) + ' Angstroms before start of spectra')
        else:
            stackstatus.append('success')

            cutstart = findval(wlen,(1030*(1+redshift)))
            cutend = findval(wlen,(1600*(1+redshift)))
            if cutend != 0 and cutstart != 0:
                cut = np.arange(cutstart,cutend)
                contfit = contfit[cut]
                flux = flux[cut]
                wlen = wlen[cut]

            wlenshift = wlen/(1+gcredshift) #+np.random.random(1)-0.5
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

            qsostacked.append(am_names[ind])
            qsostacked_z = np.append(qsostacked_z, redshift)
            qsostacked_sn = np.append(qsostacked_sn, stonforest)

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

    stackedinds.append(specstackinds)
    #add to carla stack:


    meangcwlenintpol = interpolate.interp1d(wlenhighres, meanspec, 'linear', bounds_error=False, fill_value=fillval)
    medgcwlenintpol = interpolate.interp1d(wlenhighres, medspec, 'linear', bounds_error=False, fill_value=fillval)

    if wlenhighres[0] < gcwlenmin:
        gcwlenmin = wlenhighres[0]
    if wlenhighres[-1] > gcwlenmax:
        gcwlenmax = wlenhighres[-1]

    #append stack to multidim variable:
    meanmultistore[carla, 0:] = meangcwlenintpol(wlenmultistack)
    medmultistore[carla, 0:] = medgcwlenintpol(wlenmultistack)


#average togther carla stacks:
print('total number of anti-mtach spectra = '+ str(len(qsostacked)))


gcstart = (np.abs(wlenmultistack - gcwlenmin)).argmin() +10 #tolerance for error in exact start
gcend = (np.abs(wlenmultistack - gcwlenmax)).argmin() -10
wlenmultistack = wlenmultistack[gcstart:gcend]
#cut downcols / remove empty rows
#mean
meanmultistore = meanmultistore[0:, gcstart:gcend]
meancarla = np.nanmean(meanmultistore, axis=0)
#med
medmultistore = medmultistore[0:, gcstart:gcend]
medcarla = np.nanmedian(medmultistore, axis=0)

c = 299792.458
lam = wlenmultistack
lam_em = 1215.67
vrel = c*((lam  - lam_em)/lam_em)

#ref lines from sdss
refline = open('refline.txt', 'r')
linename = []
wlenline = np.array([])
for line in refline:
    line = line.strip()
    columns = line.split()
    wlenline = np.append(wlenline, float(columns[0]))
    linename.append(columns[1])
refline.close()
linename = linename[0:12] #restrict higher WLEN lines
wlenline = wlenline[0:12]

fig = plt.figure('mean and median')
# set height ratios for sublots
gs = plt.GridSpec(2, 1, height_ratios=[1, 1])
# the fisrt subplot
ax0 = plt.subplot(gs[0])
#ax0.set_yscale("log")
line0, = ax0.plot(wlenmultistack, meancarla, color='r')
for i in range(1,12): #plot emission lines for reference
    ax0.plot(wlenline[i],meancarla[findval(wlen,wlenline[i])],'.')
ax0.set_ylim([-0.5, 1.5])
ax0.set_xlim([1100, 1300])
#the second subplot
ax1 = plt.subplot(gs[1], sharex = ax0, sharey = ax0)
line1, = ax1.plot(wlenmultistack, medcarla, color='b')
for i in range(1,12): #plot emission lines for reference
    ax1.plot(wlenline[i],medcarla[findval(wlen,wlenline[i])],'.',label = linename[i])
plt.setp(ax0.get_xticklabels(), visible=False)
# remove last tick label for the second subplot
yticks = ax1.yaxis.get_major_ticks()
yticks[-1].label1.set_visible(False)
# put leg on first subplot
ax0.legend((line0, line1), ('mean', 'median'), loc='lower left')
ax1.legend(loc='lower right')
#labels
# Set common labels
fig.text(0.5, 0.04, r'$\lambda$ ($\mathrm{\AA}$)', ha='center', va='center')
fig.text(0.06, 0.5, r'$F$ $(10^{-17}$ ergs $s^{-1}cm^{-2}\mathrm{\AA}^{-1})$', ha='center', va='center', rotation='vertical')
# remove vertical gap between subplots
plt.subplots_adjust(hspace=.0)
plt.show()

if saveoutput == True:
    wlencol = fits.Column(name='Wavelength', array = wlenmultistack, format='F')
    vrelcol = fits.Column(name='V_rel', array = vrel, format='F')
    medfluxcol = fits.Column(name='Median_Flux', array = medcarla, format='F')
    meanfluxcol = fits.Column(name='Mean_Flux', array = meancarla, format='F')
    stacksavedata = fits.BinTableHDU.from_columns([wlencol, vrelcol, medfluxcol, meanfluxcol])
    # #metadata:
    # #ext2 = qso in the stack
    qsostackedcol = fits.Column(name='QSO', array = qsostacked, format='25A')
    qso_zcol = fits.Column(name='Redshift', array = qsostacked_z, format='F')
    qso_sncol = fits.Column(name='S_N', array = qsostacked_sn, format='F')
    qsostackedsavedata = fits.BinTableHDU.from_columns([qsostackedcol, qso_zcol, qso_sncol])


    primary = fits.PrimaryHDU()
    hdul = fits.HDUList([primary, stacksavedata, qsostackedsavedata])

    outname = 'stacking/figures/Stacking data/' + runsavename+ '(control).fits'
    hdul.writeto(outname, overwrite = True)







#
