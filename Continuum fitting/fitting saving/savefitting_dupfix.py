import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate, signal
import os, time
#import fitting
import sys
#sys.path.append('C:/Users/jason/GIT/4th_year_project_git/Continuum Fitting')
sys.path.append('C:/Users/alexw/Documents/GitHub/4th_year_project_git/Continuum fitting')


#definitions for use in fitting
def findval(array,val): # find index of nearst value in array
    array = np.asarray(array)
    ind = np.abs(array - val).argmin()
    return ind

def findmax(array): #find max index
    array = np.asarray(array)
    ind = (np.abs(array - np.max(array))).argmin()
    return ind

def findmed(array): #find median index
    array = np.asarray(array)
    ind = (np.abs(array - np.median(array))).argmin()
    return ind

def findpctmax(array,pct): #find 50th percentile of top pct % in array
    ind = np.argsort(array)
    sortedarray = array[ind]
    start = int(len(sortedarray)/2 - pct/2)
    end = int(len(sortedarray)/2 + pct/2)
    selectedvals = sortedarray[start:end]
    #selectedinds = ind[-pct:] #indices of the top pct in array

    selectedinds = ind[int(-pct*0.5)]
    return selectedinds

def findpctmean(array,pct): #find median pct of data in array
#e.g if pct = 10% of array it takes 45-55 percentile
    sortedarray = np.argsort(array)
    start = int(len(sortedarray)/2 - pct/2)
    end = int(len(sortedarray)/2 + pct/2)
    selectedvals = sortedarray[start:end]
    return selectedvals



#missing files
specnames =['spec-0410-51816-0106.fits', 'spec-0410-51816-0559.fits', 'spec-0411-51817-0290.fits', 'spec-0435-51882-0630.fits', 'spec-0435-51882-0637.fits', 'spec-0516-52017-0208.fits', 'spec-0543-52017-0022.fits', 'spec-0543-52017-0152.fits', 'spec-0543-52017-0237.fits', 'spec-0543-52017-0270.fits', 'spec-0544-52201-0331.fits', 'spec-0545-52202-0464.fits', 'spec-0708-52175-0554.fits', 'spec-0834-52316-0287.fits', 'spec-0957-52398-0306.fits', 'spec-1067-52616-0268.fits', 'spec-1089-52913-0188.fits', 'spec-1159-52669-0015.fits', 'spec-1159-52669-0044.fits', 'spec-1324-53088-0620.fits', 'spec-1325-52762-0374.fits', 'spec-1512-53742-0152.fits', 'spec-1571-53174-0634.fits', 'spec-1950-53436-0452.fits', 'spec-1976-53401-0102.fits', 'spec-1976-53401-0103.fits', 'spec-2531-54572-0238.fits', 'spec-2561-54597-0630.fits.1', 'spec-2947-54533-0417.fits', 'spec-3680-55210-0128.fits.1', 'spec-3680-55210-0168.fits.1', 'spec-3680-55210-0174.fits.1', 'spec-3680-55210-0312.fits.1', 'spec-3680-55210-0314.fits.1', 'spec-3791-55501-0238.fits.1', 'spec-3791-55501-0270.fits.1', 'spec-3791-55501-0276.fits.1', 'spec-3791-55501-0310.fits.1', 'spec-3791-55501-0360.fits.1', 'spec-3791-55501-0662.fits.1', 'spec-3791-55501-0682.fits.1', 'spec-3791-55501-0704.fits.1', 'spec-3791-55501-0726.fits.1', 'spec-3791-55501-0736.fits.1', 'spec-3791-55501-0738.fits.1', 'spec-3791-55501-0792.fits.1', 'spec-3791-55501-0792.fits.2', 'spec-3791-55501-0798.fits.1', 'spec-3791-55501-0802.fits.1', 'spec-3800-55486-0362.fits.1', 'spec-3800-55486-0462.fits.1', 'spec-3800-55486-0475.fits.1', 'spec-3802-55528-0710.fits.1', 'spec-3803-55513-0470.fits.1', 'spec-3803-55513-0482.fits.1', 'spec-3803-55513-0488.fits.1', 'spec-3843-55278-0870.fits', 'spec-4457-55858-0242.fits.1', 'spec-4457-55858-0288.fits.1', 'spec-4457-55858-0288.fits.2', 'spec-4457-55858-0297.fits.1', 'spec-4457-55858-0330.fits.1', 'spec-4457-55858-0364.fits.1', 'spec-4457-55858-0364.fits.2', 'spec-4457-55858-0368.fits.1', 'spec-4457-55858-0404.fits.1', 'spec-4457-55858-0442.fits.1', 'spec-4461-55888-0030.fits.1', 'spec-4461-55888-0114.fits.1', 'spec-4461-55888-0928.fits.1', 'spec-4461-55888-0929.fits.1', 'spec-4461-55888-0978.fits.1', 'spec-4461-55888-0992.fits.1', 'spec-4461-55888-0998.fits.1', 'spec-4463-55868-0148.fits.1', 'spec-4463-55868-0185.fits.1', 'spec-4463-55868-0190.fits.1', 'spec-4463-55868-0194.fits.1', 'spec-4463-55868-0396.fits.1', 'spec-4463-55868-0398.fits.1', 'spec-4463-55868-0562.fits.1', 'spec-4463-55868-0630.fits.1', 'spec-4463-55868-0630.fits.2', 'spec-4463-55868-0640.fits.1', 'spec-4463-55868-0662.fits.1', 'spec-4463-55868-0662.fits.2', 'spec-4463-55868-0666.fits.1', 'spec-4463-55868-0666.fits.2', 'spec-4463-55868-0672.fits.1', 'spec-4463-55868-0674.fits.1', 'spec-4463-55868-0674.fits.2', 'spec-4463-55868-0676.fits.1', 'spec-4463-55868-0676.fits.2', 'spec-4463-55868-0700.fits.1', 'spec-4463-55868-0700.fits.2', 'spec-4463-55868-0702.fits.1', 'spec-4463-55868-0704.fits.1', 'spec-4463-55868-0704.fits.2', 'spec-4463-55868-0796.fits.1', 'spec-4463-55868-0796.fits.2', 'spec-4463-55868-0832.fits.1', 'spec-4463-55868-0832.fits.2', 'spec-4463-55868-0858.fits.1', 'spec-4463-55868-0860.fits.1', 'spec-4463-55868-0860.fits.2', 'spec-4463-55868-0872.fits.1', 'spec-4463-55868-0924.fits.1', 'spec-4463-55868-0924.fits.2', 'spec-4466-55857-0254.fits', 'spec-4467-55894-0538.fits.1', 'spec-4574-55621-0720.fits', 'spec-4615-55618-0838.fits.1', 'spec-4648-55673-0176.fits.1', 'spec-4648-55673-0180.fits.1', 'spec-4648-55673-0306.fits.1', 'spec-4653-55622-0492.fits.1', 'spec-4697-55660-0024.fits.1', 'spec-4697-55660-0030.fits.1', 'spec-4697-55660-0977.fits.1', 'spec-4700-55709-0278.fits.1', 'spec-4700-55709-0610.fits.1', 'spec-4700-55709-0612.fits.1', 'spec-4700-55709-0686.fits.1', 'spec-4700-55709-0688.fits.1', 'spec-4700-55709-0726.fits.1', 'spec-4700-55709-0732.fits.1', 'spec-4732-55648-0952.fits.1', 'spec-4768-55944-0510.fits.1', 'spec-4769-55931-0714.fits.1', 'spec-4769-55931-0762.fits.1', 'spec-4769-55931-0798.fits.1', 'spec-4769-55931-0848.fits.1', 'spec-4769-55931-0868.fits.1', 'spec-4769-55931-0872.fits.1', 'spec-4769-55931-0880.fits.1', 'spec-4769-55931-0894.fits.1', 'spec-4769-55931-0932.fits.1', 'spec-4769-55931-0934.fits.1', 'spec-4769-55931-0936.fits.1', 'spec-4769-55931-0938.fits.1', 'spec-4769-55931-0972.fits.1', 'spec-4849-55945-0232.fits.1', 'spec-4849-55945-0282.fits.1', 'spec-4849-55945-0303.fits.1', 'spec-4849-55945-0340.fits.1', 'spec-4849-55945-0438.fits.1', 'spec-4849-55945-0490.fits.1', 'spec-4849-55945-0640.fits.1', 'spec-4849-55945-0712.fits.1', 'spec-4850-55929-0004.fits.1', 'spec-5148-56220-0278.fits.1', 'spec-5148-56220-0312.fits.1', 'spec-5148-56220-0320.fits.1', 'spec-5148-56220-0386.fits.1', 'spec-5148-56220-0620.fits.1', 'spec-5148-56220-0690.fits.1', 'spec-5733-56575-0522.fits', 'spec-5944-56220-0122.fits.1', 'spec-5944-56220-0284.fits.1', 'spec-5983-56310-0742.fits.1', 'spec-7449-56740-0434.fits', 'spec-7876-57002-0036.fits', 'spec-8169-57071-0988.fits.1', 'spec-8193-57333-0280.fits.1', 'spec-8193-57333-0292.fits.1', 'spec-8193-57333-0306.fits.1', 'spec-8193-57333-0343.fits.1', 'spec-8193-57333-0346.fits.1', 'spec-8193-57333-0663.fits.1', 'spec-8193-57333-0740.fits.1', 'spec-8291-57391-0084.fits.1', 'spec-8291-57391-0123.fits.1', 'spec-8291-57391-0127.fits.1', 'spec-8291-57391-0128.fits.1', 'spec-8291-57391-0129.fits.1', 'spec-8291-57391-0184.fits.1', 'spec-8840-57425-0507.fits.1', 'spec-8840-57425-0561.fits.1', 'spec-8841-57425-0075.fits.1', 'spec-8841-57425-0848.fits.1', 'spec-8841-57425-0852.fits.1', 'spec-8841-57425-0871.fits.1', 'spec-8873-57427-0567.fits.1', 'spec-8874-57426-0994.fits.1']

errorcheckspec = [] #stores spec file name for those with caught forest error
dupcheckspec = []
dupcheckinstance = []

for spec in specnames: #get spec and instance it appears in pos tabel to append correct gc name
    #dupchecking:
    if spec[-1] == 's':
        errorcheckspec.append(spec)
    if spec[-1] == '1':
        dupcheckspec.append(spec[:-2])
        dupcheckinstance.append(0)
        dupcheckspec.append(spec[:-2])
        dupcheckinstance.append(1)
    if spec[-1] == '2':
        dupcheckspec.append(spec[:-2])
        dupcheckinstance.append(2)


print(len(dupcheckspec))
print(len(dupcheckinstance))




spec = 'spec-0343-51692-0145.fits'
instance = 0
#---------------------------data extraction-----------------------------#
specdirectory = 'Spectra/'+spec
#print(specdirectory)

data = fits.getdata(specdirectory,ext=1)#Read in fits spectrum data
speclen = len(data)
fitdata = fits.getdata(specdirectory,ext=2)#read in fits meta data
metasize = len(fitdata[0])
#predefine data variables
flux = np.zeros(speclen)
wlen = np.zeros(speclen)
model = np.zeros(speclen)
ivar = np.zeros(speclen)

#extract metadata
for i in range(0,speclen):
    flux[i] = data[i][0]
    ivar[i] = data[i][2]
    if ivar[i] == 0:
        ivar[i] = 0.0000001
    wlen[i] = 10**(data[i][1])
    model[i] = data[i][7]

#extract redshift of qso
if metasize == 126:
     redshift = fitdata[0][63] #accounts for diff location for sdss/eboss
else:
     redshift = fitdata[0][38]

#extract qso-gc info from positions table -> asscoaited gc and angular seperation etc...
posdata = fits.getdata('PositionsTable.fits',ext=1)
poslen = len(posdata)
clusterredshift = 0 #metadata
clusternames = 'null' #metadata
clusterseperation = 0 #metadata


findnum = 0
for j in range(0,poslen):
    #get assocaited id to match spec in folder to spec in postions table
    plate = str(posdata[j][4]).zfill(4)
    mjd = str(posdata[j][5]).zfill(5)
    fiberid = str(posdata[j][6]).zfill(4)
    specfilename = 'spec-'+plate+'-'+mjd+'-'+fiberid+'.fits'
    if specfilename == spec:
        if findnum == instance:
            clusterredshift = posdata[j][109]
            clusternames = str(posdata[j][105])
            clusterseperation = posdata[j][114]
            print(findnum)
        findnum = findnum + 1

#lymanalpha calculation:
lyalpha = 1215.67*(1+redshift)#calc lya using redshift #metadata
lyalphaind = findval(wlen,lyalpha) #metadata
gclyalpha = 1215.67*(1+clusterredshift) #metadata
gclyalphaind = findval(wlen,gclyalpha) #metdata


#s/n checking
std = (1/ivar)**0.5
stonall = np.median(flux/std) #metadata

if lyalphaind == 0: #catches error if no forest
    stonforest = 0
else:
    stonforest = np.median(flux[0:lyalphaind]/std[0:lyalphaind]) #metadata

#--------------------Continuum fitting-----------------------#

#split the spec in two about lyalpha peak for 2 fitting regions
pw = 30
intervalwlen = np.array([])
winpeak = np.array([])
forestwinnum, forestpct = 8, 0.2 #forest number must be even #metadata
otherwinnum, otherpct =  50, 0.2 #metadata

#fitv9 method
if lyalpha - wlen[0] <= 0: #if no forest just fit the other part of the spectrum
    stackstatus = 'FORESTERROR' #metadata
    step = lyalphaind
    while step <= speclen:
        window = int((speclen-lyalphaind)/otherwinnum)
        percentage = otherpct
        pct = int(window*percentage)
        windata = flux[step:(step+window)]
        winpeakind = step + findpctmean(windata,pct)
        winpeak = np.append(winpeak,flux[winpeakind])
        intervalwlen = np.append(intervalwlen,wlen[winpeakind])
        step = step + window
else: #for good spec fit accordingly
    stackstatus = 'SUCCESS' #metadata
    step = 0
    while step <= speclen:
        if step <= lyalphaind :
            winnum = forestwinnum
            window = int(lyalphaind/winnum)
            # if window < lyalphaind:
            #     window = int((speclen-lyalphaind)/otherwinnum)
            percentage = forestpct
            pct = int(window*percentage)
            windata = flux[step:(step+window)]
            winpeakind = step + findpctmean(windata,pct)#change
            winpeak = np.append(winpeak,flux[winpeakind])
            intervalwlen = np.append(intervalwlen,wlen[winpeakind])
            step = step + window
        else:
            winnum = otherwinnum
            window = int((speclen-lyalphaind)/winnum)
            percentage = otherpct
            pct = int(window*percentage)
            windata = flux[step:(step+window)]
            winpeakind = step + findpctmean(windata,pct)
            winpeak = np.append(winpeak,flux[winpeakind])
            intervalwlen = np.append(intervalwlen,wlen[winpeakind])
            step = step + window

    #pad interval with start/end value to allign correctly
intervalwlen[0] = wlen[0]
intervalwlen[-1] = wlen[-1]
intpol = interpolate.interp1d(intervalwlen, winpeak, kind=1)
contfit = intpol(wlen)
#smooth fit
if stackstatus == 'SUCCESS':
    smoothwin = 2*int(lyalphaind/forestwinnum)+1 #ensures smooth window is odd number
    if smoothwin < 10:
        smoothwin = 2*int((speclen-lyalphaind)/otherwinnum)+1
else:
    smoothwin = 2*int((speclen-lyalphaind)/otherwinnum)+1

contfit = signal.savgol_filter(contfit, smoothwin,3)
normspec = flux-contfit

wlencol = fits.Column(name='Wavelength', array=wlen, format='F')
normspeccol = fits.Column(name='Normalised_Flux', array=normspec, format='F')

specfitdata = fits.BinTableHDU.from_columns([wlencol, normspeccol])

redshiftcol = fits.Column(name='QSO_Z', array=np.array([redshift]), format='F')
stonallcol = fits.Column(name='STON_ALL', array=np.array([stonall]), format='F')
stonforestcol = fits.Column(name='STON_FOREST', array=np.array([stonforest]), format='F')
lyalphacol = fits.Column(name='QSO_LYa', array=np.array([lyalpha]), format='F')
lyalphaindcol = fits.Column(name='QSO_LYa_INDEX', array=np.array([lyalphaind]), format='K')
gcnamecol = fits.Column(name='GC_NAME', array=np.array([clusternames]), format='20A')
gcredshiftcol = fits.Column(name='GC_Z', array=np.array([clusterredshift]), format='F')
gc_qso_sepcol = fits.Column(name='GC_SEPERATION', array=np.array([clusterseperation]), format='F')
gclyalphacol = fits.Column(name='GC_LYa', array=np.array([gclyalpha]), format='F')
gclyalphaindcol = fits.Column(name='GC_LYa_INDEX', array=np.array([gclyalphaind]), format='K')
stackmsgcol = fits.Column(name='STACK_STATUS', array=np.array([stackstatus]), format='20A')

metadata = fits.BinTableHDU.from_columns([redshiftcol, stonallcol, stonforestcol, lyalphacol,
lyalphaindcol, gcnamecol, gcredshiftcol, gc_qso_sepcol, gclyalphacol, gclyalphaindcol, stackmsgcol])

primary = fits.PrimaryHDU()
hdul = fits.HDUList([primary, specfitdata, metadata])

hdul.writeto('Fitted Spectra/' + spec[0:20] +'('+clusternames+')-prefitted.fits',overwrite = True)
print('done '+spec)








    #
