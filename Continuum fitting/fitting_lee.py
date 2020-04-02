import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits
from scipy import interpolate
import os


folderpath = 'E:/spectralyalpha/BOSSLyaDR9_spectra/BOSSLyaDR9_spectra/'

#specmatch = ['spec-0343-51692-0145.fits','spec-0435-51882-0637.fits','spec-2947-54533-0417.fits','spec-3970-55591-0148.fits']

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


#given specfile name - get corresponding file in lyalpa folder on E:
#specmatch = ['spec-0343-51692-0145.fits','spec-0435-51882-0637.fits','spec-2947-54533-0417.fits','spec-3970-55591-0148.fits']
specmatch = ['spec-0343-51692-0145.fits', 'spec-0344-51693-0159.fits', 'spec-0410-51816-0106.fits', 'spec-0410-51816-0117.fits', 'spec-0410-51816-0559.fits', 'spec-0411-51817-0290.fits', 'spec-0411-51817-0399.fits', 'spec-0435-51882-0055.fits', 'spec-0435-51882-0630.fits', 'spec-0435-51882-0637.fits', 'spec-0499-51988-0560.fits', 'spec-0508-52366-0520.fits', 'spec-0516-52017-0208.fits', 'spec-0543-52017-0022.fits', 'spec-0543-52017-0152.fits', 'spec-0543-52017-0237.fits', 'spec-0543-52017-0270.fits', 'spec-0544-52201-0331.fits', 'spec-0544-52201-0419.fits', 'spec-0545-52202-0464.fits', 'spec-0708-52175-0153.fits', 'spec-0708-52175-0554.fits', 'spec-0709-52205-0170.fits', 'spec-0709-52205-0483.fits', 'spec-0709-52205-0516.fits', 'spec-0724-52254-0113.fits', 'spec-0725-52258-0182.fits', 'spec-0725-52258-0295.fits', 'spec-0832-52312-0113.fits', 'spec-0832-52312-0501.fits', 'spec-0832-52312-0566.fits', 'spec-0833-52314-0084.fits', 'spec-0833-52314-0194.fits', 'spec-0833-52314-0450.fits', 'spec-0834-52316-0243.fits', 'spec-0834-52316-0287.fits', 'spec-0834-52316-0292.fits', 'spec-0834-52316-0391.fits', 'spec-0862-52325-0386.fits', 'spec-0915-52443-0111.fits', 'spec-0937-52707-0232.fits', 'spec-0957-52398-0306.fits', 'spec-1006-52708-0462.fits', 'spec-1006-52708-0497.fits', 'spec-1067-52616-0196.fits', 'spec-1067-52616-0205.fits', 'spec-1067-52616-0268.fits', 'spec-1067-52616-0292.fits', 'spec-1067-52616-0433.fits', 'spec-1068-52614-0634.fits', 'spec-1089-52913-0188.fits', 'spec-1159-52669-0015.fits', 'spec-1159-52669-0044.fits', 'spec-1203-52669-0576.fits', 'spec-1267-52932-0264.fits', 'spec-1316-52790-0354.fits', 'spec-1324-53088-0620.fits', 'spec-1325-52762-0374.fits', 'spec-1457-53116-0205.fits', 'spec-1463-53063-0029.fits', 'spec-1512-53742-0030.fits', 'spec-1512-53742-0152.fits', 'spec-1513-53741-0170.fits', 'spec-1513-53741-0216.fits', 'spec-1513-53741-0222.fits', 'spec-1513-53741-0437.fits', 'spec-1513-53741-0463.fits', 'spec-1513-53741-0466.fits', 'spec-1563-52941-0228.fits', 'spec-1567-53172-0058.fits', 'spec-1567-53172-0123.fits', 'spec-1571-53174-0634.fits', 'spec-1586-52945-0179.fits', 'spec-1732-53501-0205.fits', 'spec-1737-53055-0335.fits', 'spec-1852-53534-0401.fits', 'spec-1911-53295-0044.fits', 'spec-1911-53295-0212.fits', 'spec-1950-53436-0452.fits', 'spec-1971-53472-0578.fits', 'spec-1976-53401-0052.fits', 'spec-1976-53401-0102.fits', 'spec-1976-53401-0103.fits', 'spec-1977-53475-0221.fits', 'spec-2063-53359-0351.fits', 'spec-2068-53386-0305.fits', 'spec-2168-53886-0222.fits', 'spec-2278-53711-0290.fits', 'spec-2424-54448-0256.fits', 'spec-2443-54082-0178.fits', 'spec-2529-54585-0635.fits', 'spec-2531-54572-0233.fits', 'spec-2531-54572-0238.fits', 'spec-2531-54572-0246.fits', 'spec-2561-54597-0630.fits', 'spec-2610-54476-0341.fits', 'spec-2782-54592-0488.fits', 'spec-2782-54592-0512.fits', 'spec-2947-54533-0417.fits', 'spec-3606-55182-0033.fits', 'spec-3606-55182-0090.fits', 'spec-3606-55182-0134.fits', 'spec-3609-55201-0706.fits', 'spec-3609-55201-0796.fits', 'spec-3609-55201-0834.fits']
specmatchlya =[]
for j in range(0,len(specmatch)):
    specfilename = specmatch[j]
    out = 'null'
    #finds matching spec
    for i in range(0,matchlen):
        check = specpath[i]
        if check == specfilename:
            out = speclyapath[i]
    specmatchlya.append(out)
    print('For DR14 '+specfilename+' corresponding spec in lya continuum fit is '+out[-28:])



#show matching spec
for i in range(len(specmatchlya)):
    if specmatchlya[i] != 'null':
        data = fits.getdata(specmatchlya[i],ext=1)#import fits image
        speclen = len(data)
        flux = np.zeros(speclen)
        wlen = np.zeros(speclen)
        model = np.zeros(speclen)
        ivar = np.zeros(speclen)
        cont = np.zeros(speclen)

        for i in range(0,speclen):
            flux[i] = data[i][0]
            ivar[i] = data[i][2]
            if ivar[i] == 0:
                ivar[i] = 0.00001
            wlen[i] = 10**(data[i][1])
            model[i] = data[i][7]
            cont[i] = data[i][11]

        contrem =flux-cont

        plt.figure('flux and fitting')
        plt.plot(wlen,flux, label='flux')
        plt.plot(wlen,cont, label='cont fit')
        plt.legend()

        plt.figure('continuum subtraction')
        plt.plot(wlen,contrem, label='continuum removed')
        plt.legend()
        plt.show()












 #
