# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 14:44:03 2024

@author: cammy
"""

from astropy.io import fits
from os import path
import matplotlib.pyplot as plt
import numpy as np
import h5py


from os import path
import numpy as np
import matplotlib.pyplot as plt

import vorbin
from vorbin.voronoi_2d_binning import voronoi_2d_binning


   
from plotbin import display_pixels
from plotbin.display_pixels import display_pixels
from plotbin.display_bins import display_bins
from ppxf.ppxf import ppxf
import ppxf.ppxf_util as util
import ppxf.miles_util as lib

hf = h5py.File('/fred/oz059/cammy/data_processing/MAUVE-project/vorbin_NGC4694.h5', 'r')
hf.keys()
binNum=hf.get('binNum')
nPixels=hf.get('nPixels')
scale=hf.get('scale')
sn=hf.get('sn')
x_bar=hf.get('x_bar')
x_gen=hf.get('x_gen')
y_bar=hf.get('y_bar')
y_gen=hf.get('y_gen')

binNum=np.array(binNum)
nPixels=np.array(nPixels)
scale=np.array(scale)
sn=np.array(sn)
x_bar=np.array(x_bar)
x_gen=np.array(x_gen)
y_bar=np.array(y_bar)
y_gen=np.array(y_gen)
hf.close()


hdul = fits.open('/fred/oz059/')

lamda_ref=hdul[1].header['CRVal3']
nlamda=hdul[1].header['Naxis3']
refpix=hdul[1].header['CRpix']
deltlamba=hdul[1].header['CD3_3']

lamGal=lamda_ref+(nlamda-refpix)*dellamba

spec=hdul[0].data
variance=hdul[1].data

loc=np.where(np.nanmedian(spec,axis=0)>0.)

spectra_flat=spec[:,loc[0],loc[1]]
variance_flat=variance[:,loc[0],loc[1]]


nbins=max(binNum)+1
bin_spec=np.zeros([nlamda,nbins])
bin_var=np.empty_like(bin_spec)

for i in range(nbins):
    w, =np.where(binNum==i)
    bin_spec[:,i]=np.sum(spectra_flat[:,w])
    bin_var[:,i]=np.sum(variance_flat[:,w])
    binsize=len(w)
    
    
fig, axs = plt.subplots(2,2,figsize=(19.5,5),layout='constrained')
for i in range(2):
    inter=np.random.randint(0,nbins)
    axs[i,0].plot(lamGal,bin_spec[:,inter])
    axs[i,0].set_title(f'{x_bar[inter]},{y_bar[inter]}')
    inter=np.random.randint(0,nbins)
    axs[i,1].set_title(f'{x_bar[inter]},{y_bar[inter]}')
    axs[i,1].plot(lamGal,bin_spec[:,inter])
plt.savefig('binned_spectra_check.png')    




lamRange1 = numpy.array([numpy.min(lamGal),numpy.max(lamGal)])

# We log-rebin each row of the galaxy spectra in turn, so lets do that in
# a loop, and store the log-binned spectra in a new 2D array
for i in range(0,nbins):

  # Run the i'th spectrum through the log-rebin routine, storing in 'galaxy'
  # logLam1 will be the new wavelength limits
  # velscale will be the new pixel size, in km/s
  galaxy, logLam1, velscale = util.log_rebin(lamRange1, bin_spec[:,i])
  variance, logLam1, velscale = util.log_rebin(lamRange1, bin_var[:,i])

  # Now that we know how long the output spectra will be, we can create the
  # 2D output array for the log-binned galaxy spectra
  if i ==0:
    log_gal = numpy.zeros((galaxy.shape[0],nbins))
    log_var = numpy.zeros((variance.shape[0],nbins))

  # Now insert the rebinned spectrum into the 2D array in the correct row
  log_gal[:,i] = galaxy
  log_var[:,i] = variance



print("Log-rebinned data array shape:",log_gal.shape)
print("Velscale: ",velscale, "km/s")

fig, axs = plt.subplots(2,2,figsize=(19.5,5),layout='constrained')
for i in range(2):
    inter=np.random.randint(0,nbins)
    axs[i,0].plot(loglam1,log_gal[:,inter])
    axs[i,0].set_title(f'{x_bar[inter]},{y_bar[inter]}')
    inter=np.random.randint(0,nbins)
    axs[i,1].set_title(f'{x_bar[inter]},{y_bar[inter]}')
    axs[i,1].plot(loglam1,log_gal[:,inter])
plt.savefig('log_spectra_check.png')   






