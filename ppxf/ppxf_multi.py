#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 11:01:23 2024

@author: 48105686
"""
from astropy.io import fits

import matplotlib.pyplot as plt
import numpy as np
import h5py   
from plotbin import display_pixels
from plotbin.display_pixels import display_pixels
from plotbin.display_bins import display_bins


from ppxf.ppxf import ppxf
import ppxf.ppxf_util as util
import ppxf.sps_util as lib
from time import perf_counter as clock
from urllib import request
from multiprocessing import Pool
import os

os.environ['OPENBLAS_NUM_THREADS'] = '1'

def create_template_file():
    




    dir_path='C:/Users/cammy/OneDrive/Documents/astronomy/ppxf/E-MILES'
    dir_path='/Users/48105686/Public/Drop Box/astronomy/ppxf/E-MILES'

    res = []

# Iterate directory
    for path in os.listdir(dir_path):
    # check if current path is a file
        if os.path.isfile(os.path.join(dir_path, path)):
            res.append(path)
            check=0
            for i in range(len(res)):
                if res[i]=='spectra_e-miles_1.3.npz':
                    check=1


    from os import path
    if check != 1:
        print('creating file')
        hdul=fits.open(dir_path+'/'+res[68])
        data=hdul[0].data
        nlam=hdul[0].header['NAXIS1']
        lam_start=hdul[0].header['CRVAL1']
        lam_delt=hdul[0].header['CDELT1']
        lam=lam_start+(np.arange(nlam)*lam_delt)
        temp=np.zeros([nlam,len(res)])
        
        temp_m=np.zeros_like(res)
        
        temp_z=np.zeros_like(res)
        temp_a=np.zeros_like(res)
        
        #print(res[1],res[0][3:7],res[0][9:13],res[0][14:21])
        
        for i in range(len(res)):
            #print(float(res[i][3:7]),res[i][8],res[i][9:13],float(res[i][14:20]))
            print(i,res[i])
            temp_m[i]=float(res[i][3:7])
            if res[i][8]=='m':
                temp_z[i]=-1*float(res[i][9:13])
            else:
                temp_z[i]=float(res[i][9:13])
            temp_a[i]=float(res[i][14:21])
                #print(float(res[i][14:21]),res[i][14:21])
            hdul=fits.open(dir_path+'/'+res[i])
                
            temp[:,i]=hdul[0].data
                    
        metals=np.unique(temp_z)
        ages=np.unique(temp_a)
        print(temp_a)
        masses=np.zeros([len(ages),len(metals)])
        print(ages,metals)
        masses=np.full(masses.shape,1.3)
        print(masses)
        templates=np.zeros([nlam,len(ages),len(metals)])
        print(temp.shape,metals.size,ages.size)
        for i in range(len(metals)):
            for j in range(len(ages)):
                templates[:,j,i]
                templates[:,j,i]=temp[:,len(ages)*i+j]
                            
        FWHM=np.zeros(lam.shape)
                            
        fwhm_3060=3
        fwhm_3541=5
        fwhm_8950=2.51
        fwhm_50000=c*lam_delt/lam
                            
        for i in range(len(lam)):
            #print(i)
            if lam[i]<=3060.8:
                FWHM[i]=fwhm_3060
            elif lam[i]>3060 and lam[i]<=3541:
                FWHM[i]=fwhm_3541
            elif lam[i]>3541 and lam[i]<=8950:
                FWHM[i]=fwhm_8950
            else:
                FWHM[i]=fwhm_50000[i]
                                                
        np.savez_compressed('spectra_e-miles_1.3.npz', templates=templates, masses=masses, 
                                                                    lam=lam, ages=ages, metals=metals, fwhm=FWHM)
                                                
def get_data():
    hf = h5py.File('/Users/48105686/Public/Drop Box/astronomy/fits/vorbin_NGC4694.h5','r')
    binNum=np.array(hf.get('binNum'))
    x_bar=np.array(hf.get('x_bar'))
    y_bar=np.array(hf.get('y_bar'))
    x_gen=np.array(hf.get('x_gen'))
    y_gen=np.array(hf.get('y_gen'))
    sn=np.array(hf.get('sn'))
    nPixels=np.array(hf.get('nPixels'))
    scale=np.array(hf.get('scale'))
    target_sn=np.array(hf.get('target_sn'))
    lamGal=np.array(hf.get('lamGal'))
    bin_spec=np.array(hf.get('bin_spec'))
    bin_var=np.array(hf.get('bin_var'))
    binsize=np.array(hf.get('bin_size'))
    print(hf.keys())
    hf.close()
    return binNum,x_bar,y_bar,x_gen,y_gen,sn,nPixels,scale,target_sn,lamGal,bin_spec,bin_var,binsize



def log_rebin(lamGal,bin_spec,bin_var,nbins):
    lamRange1 = np.array([np.min(lamGal),np.max(lamGal)])

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
            log_gal = np.zeros((galaxy.shape[0],nbins))
            log_var = np.zeros((variance.shape[0],nbins))

      # Now insert the rebinned spectrum into the 2D array in the correct row
        log_gal[:,i] = galaxy
        log_var[:,i] = variance
    return log_gal,log_var,logLam1,velscale


############################################################################

def summed_ppxf(filename,velscale,fwhm_gal,lamGal,log_gal,log_var,logLam1):
    
    lam_range_temp = [lamGal[0]/1.02, lamGal[-1]*1.02]
   
    #templates=sps.templates
    #print(templates.shape)
    # Compute a mask for gas emission lines
    #goodPixels = util.determine_goodpixels(logLam1, lam_range_temp, redshift)
    print(sps.templates.shape,np.reshape(sps.age_grid,-1))
    
    stars_templates = sps.templates.reshape(sps.templates.shape[0], -1)
    ages=np.reshape(sps.age_grid,-1)
    metals=np.reshape(sps.metal_grid,-1)
    gas_templates, gas_names, line_wave = util.emission_lines(
    sps.ln_lam_temp, lam_range_gal, fwhm_gal, tie_balmer=tie_balmer,
        limit_doublets=limit_doublets)
        # Here the actual fit starts. The best fit is plotted on the screen. Gas
        # emission lines are excluded from the pPXF fit using the GOODPIXELS
        # keyword.
        
    templates = np.column_stack([stars_templates, gas_templates])
        
    print(templates,stars_templates.shape)

    n_temps = stars_templates.shape[1]
    n_forbidden = np.sum(["[" in a for a in gas_names])  # forbidden lines contain "[*]"
    n_balmer = len(gas_names) - n_forbidden
    
        # Assign component=0 to the stellar templates, component=1 to the Balmer
        # gas emission lines templates and component=2 to the forbidden lines.
    component = [0]*n_temps + [1]*n_balmer + [2]*n_forbidden
    gas_component = np.array(component) > 0  # gas_component=True for gas templates

    gas_reddening = 0 if tie_balmer else None

    # Fit (V, sig, h3, h4) moments=4 for the stars
    # and (V, sig) moments=2 for the two gas kinematic components
    moments = np.array([4, 2, 2])

    # Adopt the same starting value for the stars and the two gas components



    vel = c*np.log(1 + redshift)   # eq.(8) of Cappellari (2017, MNRAS)
    start = [vel, 200.]  # (km/s), starting guess for [V, sigma]
    start = np.array([start, start, start])
    print(len(start[1]),moments[1])
    t = clock()
    summed_spec=np.nansum(log_gal,axis=1)
    
    summed_var=np.nansum(log_var,axis=1)
    print(summed_var[3200:])
    pp = ppxf(templates, summed_spec, np.sqrt(summed_var), velscale, start, plot=True, lam=np.exp(logLam1),
          lam_temp=sps.lam_temp, degree=4,component=component,
          gas_component=gas_component, gas_names=gas_names,
          gas_reddening=gas_reddening,moments=moments)
#print(pp.gas_bestfit,pp.gas_flux)
#print(np.array(np.where(pp.weights.reshape(sps.n_ages,sps.n_metals)>0)).T)
#print(pp.weights[:pp.weights.shape[0]-len(gas_names)])

    weights=np.array(pp.weights)
    print(weights,np.where(weights>0))
#temp_temp=np.zeros([sps.lam_temp.size,sps.n_ages,sps.n_metals])
#for i in range(sps.n_ages):
#    for j in range(sps.n_metals):
#        if weights[i,j]>0:
#            temp_temp[:,i,j]=templates[:,i,j]
    temp_temp2=templates[:,np.where(weights>0)[0]]
    ages_temp2=ages[np.where(weights[:weights.size-len(gas_names)]>0)[0]]
    metals_temp2=metals[np.where(weights[:weights.size-len(gas_names)]>0)[0]]

    print('cunt',temp_temp2.shape)
    return temp_temp2,ages_temp2,metals_temp2










def ppxf_f(index):
    pp = ppxf(temp_temp2, log_gal[:,index], np.sqrt(log_var[:index]), velscale, start, plot=False, lam=np.exp(logLam1),
              lam_temp=sps.lam_temp, degree=4,component=component,
              gas_component=gas_component, gas_names=gas_names,
              gas_reddening=gas_reddening,moments=moments)
    
    V,rms,h3,h4=pp.sol[0]
    Vhb,rmshb=pp.sol[1]
    Vo3,rmso3=pp.sol[2]
    eV,erms,eh3,eh4=pp.error[0]*np.sqrt(pp.chi2)
    chi2=pp.chi2
    goodpixels=pp.goodpixels
    #print(pp.weights.reshape(25,6))
    weights=pp.weights
    redshift_fit = (1 + redshift_0)*np.exp(V/c) - 1  # eq. (5c) C22
    redshift_err = (1 + redshift_fit)*eV/c    
    bestfits=pp.bestfit
    for j in range(len(logLam1)):
   ##print(np.any(logLam1[goodpixels]==logLam1[i]))
       if np.any(logLam1[goodpixels]==logLam1[j])==True:
           masked_pixels[j]=False
    
    if index % 100 ==0 :
        print(index)
    
    return masked_pixels,V,rms,h3,h4,eV,erms,eh3,eh4,chi2,redshift_fit,redshift_err,Vhb,Vo3,rmshb,rmso3,weights,bestfits






#trunc,=np.where(log_gal[:,499]>0)#np.where((lamGal>0) & (lamGal<100000))
#trunc_lamGal=np.array(lamGal[trunc])
#trunc_loglam=np.array(logLam1[trunc])
#trunc_log_gal=np.array(log_gal[trunc,:])
#trunc_log_var=np.array(log_var[trunc,:])
redshift_0=0.00
redshift=0.003869
fwhm_gal=2.65
filename='spectra_e-miles_1.3.npz'
tie_balmer=False
limit_doublets=False


c = 299792.458



create_template_file()
print(5)
binNum,x_bar,y_bar,x_gen,y_gen,sn,nPixels,scale,target_sn,lamGal,bin_spec,bin_var,binsize=get_data()
nbins=np.max(binNum)+1






log_gal,log_var,logLam1,velscale= log_rebin(lamGal,bin_spec,bin_var,nbins)

lam_range_temp = [lamGal[0]/1.02, lamGal[-1]*1.02]
lam_range_gal = np.array([np.min(lamGal), np.max(lamGal)])/(1 + redshift)

sps = lib.sps_lib(filename, velscale, fwhm_gal, wave_range=lam_range_temp)


temp_temp2,ages_temp2,metals_temp2=summed_ppxf(filename,velscale,fwhm_gal,lamGal,log_gal,log_var,logLam1)


gas_templates, gas_names, line_wave = util.emission_lines(
sps.ln_lam_temp, lam_range_gal, fwhm_gal, tie_balmer=tie_balmer,
    limit_doublets=limit_doublets)

n_temps = temp_temp2.shape[1]-len(gas_names)
n_forbidden = np.sum(["[" in a for a in gas_names])  # forbidden lines contain "[*]"
n_balmer = len(gas_names) - n_forbidden

# Assign component=0 to the stellar templates, component=1 to the Balmer
# gas emission lines templates and component=2 to the forbidden lines.
component = [0]*n_temps + [1]*n_balmer + [2]*n_forbidden
gas_component = np.array(component) > 0  # gas_component=True for gas templates

gas_reddening = 0 if tie_balmer else None



vel = c*np.log(1 + redshift)   # eq.(8) of Cappellari (2017, MNRAS)
start = [vel, 200.]  # (km/s), starting guess for [V, sigma]
start = np.array([start, start, start])
items=np.zeros([nbins,log_gal.shape[0],log_gal.shape[0]])
moments = np.array([4, 2, 2])

masked_pixels=np.full(log_gal.shape,True)
#plt.plot(logLam1,log_var[:,x])
bestfits=np.zeros(log_gal.shape)
V,rms,h3,h4,eV,erms,eh3,eh4,chi2,redshift_fit,redshift_err = np.zeros(nbins),np.zeros(nbins),np.zeros(nbins),np.zeros(nbins),np.zeros(nbins),np.zeros(nbins),np.zeros(nbins),np.zeros(nbins),np.zeros(nbins),np.zeros(nbins),np.zeros(nbins)
Vhb,Vo3,rmshb,rmso3=np.zeros(nbins),np.zeros(nbins),np.zeros(nbins),np.zeros(nbins)
weights= np.zeros([nbins,temp_temp2[0].shape[0]])




for i in range(nbins):
    items[i,0]=log_gal[:,i]
    items[i,1]=log_var[:,i]
        
items=np.arange(nbins)
i=0
with Pool() as pool:
    for result in pool.map(ppxf_f,items):
        masked_pixels[:,i],V[i],rms[i],h3[i],h4[i],eV[i],erms[i],eh3[i],eh4[i],chi2[i],redshift_fit[i],redshift_err[i],Vhb[i],Vo3[i],rmshb[i],rmso3[i],weights[i],bestfits[i]=result
        i+=1
 

print(V)
f = h5py.File('/Users/48105686/Public/Drop Box/astronomy/ppxf/ppxf_out_NGC4694.h5','w')

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
    
f.create_dataset("redshift", data=redshift_fit)
f.create_dataset("eredshift", data=redshift_err)
    
f.create_dataset("chi2", data=chi2)
f.create_dataset("spectra", data=log_gal)
    
f.create_dataset("variance", data=log_var)
f.create_dataset("bestfits", data=bestfits)
    
f.create_dataset("binNum", data=binNum)
f.create_dataset("masked_pixels", data=masked_pixels)
    
f.create_dataset("weights", data=weights)
f.create_dataset("lamGal", data=lamGal)
f.create_dataset("logLam1", data=logLam1)
f.create_dataset('ages',data=ages_temp2)
f.create_dataset('metals',data=metals_temp2)
    
f.close()    