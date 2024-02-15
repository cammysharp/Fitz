# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 16:01:56 2024

@author: cammy
"""
import h5py
import numpy as np

gal_name='NGC4694'
hf = h5py.File(f'ppxf_gas_out_{gal_name}_{1}.h5','r')
   
V=np.array(hf.get('V'))
Vhb=np.array(hf.get('Vhb'))
Vo3=np.array(hf.get('Vo3'))
bestfits=np.array(hf.get('bestfits'))
binNum=np.array(hf.get('binNum'))
chi2=np.array(hf.get('chi2'))
eV=np.array(hf.get('eV'))
eh3=np.array(hf.get('eh3'))
eh4=np.array(hf.get('eh4'))
eredshift=np.array(hf.get('eredshift'))
erms=np.array(hf.get('erms'))
h3=np.array(hf.get('h3'))
h4=np.array(hf.get('h4'))
lamGal=np.array(hf.get('lamGal'))
logLam1=np.array(hf.get('logLam1'))
masked_pixels=np.array(hf.get('masked_pixels'))
redshift=np.array(hf.get('redshift'))
rms=np.array(hf.get('rms'))
rmshb=np.array(hf.get('rmshb'))
rmso3=np.array(hf.get('rmso3'))
spectra=np.array(hf.get('spectra'))
variance=np.array(hf.get('variance'))
weights=np.array(hf.get('weights'))
x_bar=np.array(hf.get('x_bar'))
y_bar=np.array(hf.get('y_bar'))

max_rank=8
nbins=max(binNum)+1
#max_rank_range=[(max_rank-1)*nbins/max_rank,(max_rank)*nbins/max_rank]


for i in range(max_rank):
    rank_range=[(i-1)*nbins/max_rank,(i)*nbins/max_rank]
    hf = h5py.File(f'ppxf_gas_out_{gal_name}_{i}.h5','r')
    
    V[rank_range[0]:rank_range[1]]=np.array(hf.get('V'))[rank_range[0]:rank_range[1]]
    Vhb[rank_range[0]:rank_range[1]]=np.array(hf.get('Vhb'))[rank_range[0]:rank_range[1]]
    Vo3[rank_range[0]:rank_range[1]]=np.array(hf.get('Vo3'))[rank_range[0]:rank_range[1]]
    x_bar[rank_range[0]:rank_range[1]]=np.array(hf.get('x_bar'))[rank_range[0]:rank_range[1]]
    y_bar[rank_range[0]:rank_range[1]]=np.array(hf.get('y_bar'))[rank_range[0]:rank_range[1]]
    chi2[rank_range[0]:rank_range[1]]=np.array(hf.get('chi2'))[rank_range[0]:rank_range[1]]
    eV[rank_range[0]:rank_range[1]]=np.array(hf.get('eV'))[rank_range[0]:rank_range[1]]
    eh3[rank_range[0]:rank_range[1]]=np.array(hf.get('eh3'))[rank_range[0]:rank_range[1]]
    eh4[rank_range[0]:rank_range[1]]=np.array(hf.get('eh4'))[rank_range[0]:rank_range[1]]
    eredshift[rank_range[0]:rank_range[1]]=np.array(hf.get('eredshift'))[rank_range[0]:rank_range[1]]
    erms[rank_range[0]:rank_range[1]]=np.array(hf.get('erms'))[rank_range[0]:rank_range[1]]
    h3[rank_range[0]:rank_range[1]]=np.array(hf.get('h3'))[rank_range[0]:rank_range[1]]
    h4[rank_range[0]:rank_range[1]]=np.array(hf.get('h4'))[rank_range[0]:rank_range[1]]
    redshift[rank_range[0]:rank_range[1]]=np.array(hf.get('redshift'))[rank_range[0]:rank_range[1]]
    rms[rank_range[0]:rank_range[1]]=np.array(hf.get('rms'))[rank_range[0]:rank_range[1]]
    rmshb[rank_range[0]:rank_range[1]]=np.array(hf.get('rmshb'))[rank_range[0]:rank_range[1]]
    rmso3[rank_range[0]:rank_range[1]]=np.array(hf.get('rmso3'))[rank_range[0]:rank_range[1]]
    
    weights[:,rank_range[0]:rank_range[1]]=np.array(hf.get('weights'))[rank_range[0]:rank_range[1]]
    
    masked_pixels[:,rank_range[0]:rank_range[1]]=np.array(hf.get('masked_pixels'))[rank_range[0]:rank_range[1]]
    bestfits[:,rank_range[0]:rank_range[1]]=np.array(hf.get('bestfits'))[rank_range[0]:rank_range[1]]
    spectra[:,rank_range[0]:rank_range[1]]=np.array(hf.get('spectra'))[rank_range[0]:rank_range[1]]
    variance[:,rank_range[0]:rank_range[1]]=np.array(hf.get('variance'))[rank_range[0]:rank_range[1]]
    
    hf.close()

f = h5py.File(f'C:/Users/cammy/OneDrive/Documents/astronomy/ppxf/ppxf_out_{gal_name}.h5','w')

f.create_dataset("x_bar", data=x_bar)
f.create_dataset("y_bar", data=y_bar)


f.create_dataset("Vhb", data=Vhb)
f.create_dataset("rmshb", data=rmshb)

f.create_dataset("Vo3", data=Vo3)
f.create_dataset("rmso3", data=rmso3)

f.create_dataset("V", data=V)
f.create_dataset("rms", data=rms)

f.create_dataset("h3", data=h3)
f.create_dataset("h4", data=h4)

f.create_dataset("eV", data=eV)
f.create_dataset("erms", data=erms)

f.create_dataset("eh3", data=eh3)
f.create_dataset("eh4", data=eh4)

f.create_dataset("redshift", data=redshift)
f.create_dataset("eredshift", data=eredshift)

f.create_dataset("chi2", data=chi2)
f.create_dataset("spectra", data=spectra)

f.create_dataset("variance", data=variance)
f.create_dataset("bestfits", data=bestfits)

f.create_dataset("binNum", data=binNum)
f.create_dataset("masked_pixels", data=masked_pixels)

f.create_dataset("weights", data=weights)
f.create_dataset("lamGal", data=lamGal)
f.create_dataset("logLam1", data=logLam1)


hf.close()





