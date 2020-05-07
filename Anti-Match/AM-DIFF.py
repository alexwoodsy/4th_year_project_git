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

amoldtabledata = fits.getdata('Anti-Match/Anti-Match-Metadata.fits',ext=1)
amoldtablelen = len(amoldtabledata)

metaqsoname = amoldtabledata.field(1)
metara = amoldtabledata.field(2)
metadec = amoldtabledata.field(3)
metaplate = amoldtabledata.field(5).astype(int)
metamjd = amoldtabledata.field(6).astype(int)
metafiberid = amoldtabledata.field(7).astype(int)
metaz = amoldtabledata.field(9)


#imports the spectra anitmatched
specnames = next(os.walk('E:/AM-Prefitted Spectra'))[2]
spectot = len(specnames)
print(spectot)

metaind = np.array([]).astype(int)


for spec in specnames:
    count = 0
    for i in range(0,amoldtablelen):
        oldname = 'AM-spec-'+str(metaplate[i]).zfill(4)+'-'+str(metamjd[i]).zfill(5)+'-'+str(metafiberid[i]).zfill(4)+'-Preffited.fits'
        if spec == oldname and count == 0:
            print('matched '+ spec)
            metaind = np.append(metaind,i)
            count = count + 1

print(len(metaind))



#select non dups:
metaqsoname = metaqsoname[metaind]
metara = metara[metaind]
metadec = metadec[metaind]
metaplate = metaplate[metaind]
metamjd = metamjd[metaind]
metafiberid = metafiberid[metaind]
metaz = metaz[metaind]



#put in from_columns
qsonamecol = fits.Column(name='QSO_NAME', array = metaqsoname, format='20A')
racol = fits.Column(name='RA', array = metara, format='F')
deccol = fits.Column(name='DEC', array = metadec, format='F')
platecol = fits.Column(name='PLATE', array = metaplate, format='K')
mjdcol = fits.Column(name='MJD', array = metamjd, format='K')
fiberidcol = fits.Column(name='FIBERID', array = metafiberid, format='K')
zcol = fits.Column(name='Z', array = metaz, format='F')




metadata = fits.BinTableHDU.from_columns([qsonamecol ,racol ,deccol ,platecol ,mjdcol ,fiberidcol ,zcol])

primary = fits.PrimaryHDU()
hdul = fits.HDUList([primary, metadata])

hdul.writeto('Anti-Match/AM-nodups.fits',overwrite = True)
