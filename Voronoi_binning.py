# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 14:24:46 2024

@author: cammy
"""

from astropy.io import fits
from os import path
import matplotlib.pyplot as plt
import numpy as np
import h5py
import vorbin
from vorbin.voronoi_2d_binning import voronoi_2d_binning

hf = h5py.File('s2n_map.h5', 'r')

n1=hf.get('sn')
n1=np.asarray(n1)



x=np.zeros(n1[0,:,:].size)
y=np.zeros(n1[0,:,:].size)

for i in range(l):
    y[i*m:i*m+m]=l-i



for i in range(m):
    for j in range(l):
        x[j+i*l]=j
       
signal=np.zeros(n1[0,:,:].size)   
noise=np.zeros(n1[1,:,:].size)   

for i in range(m):
    signal[i*l:i*l+l]=n1[0,i,::-1]
    noise[i*l:i*l+l]=n1[1,i,::-1]
x=x[::-1]

loc=np.where(noise==0)
noise[loc]=noise[loc]+1

target_sn = 50.0

binNum, x_gen, y_gen, x_bar, y_bar, sn, nPixels, scale = voronoi_2d_binning(
        x, y, signal, noise, target_sn, plot=1, quiet=1,pixelsize=0.025)


hf = h5py.File('vorbin_NGC4694.h5','r')
hf.create_dataset('binNum',data=binNum)
hf.create_dataset('x_gen',data=x_gen)
hf.create_dataset('y_gen',data=y_gen)
hf.create_dataset('x_bar',data=x_bar)
hf.create_dataset('y_bar',data=y_bar)
hf.create_dataset('sn',data=sn)
hf.create_dataset('nPixels',data=nPixels)
hf.create_dataset('scale',data=scale)
hf.close

