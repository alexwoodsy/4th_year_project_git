#metadata extract
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.io import fits

fitdata = fits.getdata('Spectra/spec-7820-56984-0276.fits',ext=2)#import fits image

print(len(fitdata[0]))
