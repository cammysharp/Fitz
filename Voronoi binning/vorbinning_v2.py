# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 15:14:04 2024

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

def transform(x,y,angle):
    theta=np.radians(angle)
    xm=x*np.cos(theta)-y*np.sin(theta)
    ym=x*np.sin(theta)+y*np.cos(theta)
    return xm,ym


file_dir = path.dirname(path.realpath('sn_map_NGC4694.h5'))


hf = h5py.File(file_dir+'\\sn_map_NGC4694.h5', 'r')
centre=hf.get('centre')
n1=hf.get('sn')
rotation=hf.get('rot')
hf.close
n1=np.array(n1)
centre=np.array(centre)
yc,xc=centre
rotation=float(rotation)


X,Y=np.meshgrid(0.2*np.linspace(0-xc,l-xc,l),0.2*np.linspace(0-yc,m-yc,m))
X,Y=transform(X,Y,rotation)
mask=np.where(n1[0]>0)
binNum, x_gen, y_gen, x_bar, y_bar, sn, nPixels, scale = voronoi_2d_binning(
    X[mask], Y[mask], n1[0][mask], n1[1][mask], target_sn, plot=1, quiet=0,pixelsize=0.2)

plt.tight_layout()
plt.pause(1)
plt.savefig('Vorbin_plots_NGC4694.png')

hf= h5py.File('vorbin_NGC4694.h5','w')
hf.create_dataset('binNum',data=binNum)
hf.create_dataset('x_gen',data=x_gen)
hf.create_dataset('y_gen',data=y_gen)
hf.create_dataset('x_bar',data=x_bar)
hf.create_dataset('y_bar',data=y_bar)
hf.create_dataset('sn',data=sn)
hf.create_dataset('nPixels',data=nPixels)
hf.create_dataset('scale',data=scale)
hf.close()
