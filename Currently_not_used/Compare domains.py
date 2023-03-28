#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 18:10:09 2023

@author: tijmen
"""
import os
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np

os.chdir('/')

import Thesis_Functions.calculations as calculations
import Thesis_Functions.data as data
import Thesis_Functions.plotting as plotting


fig_save_directory = '/Users/tijmen/Desktop/Figures_Thesis/'
remapdir = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Data/Remap/'


RC4_GR = xr.open_dataset(remapdir+'Greenland_EASE.nc')
RC4_ARC = xr.open_dataset(remapdir+'ARCTIC_EASE.nc')

RC4_GR = RC4_GR['botflx'].isel(time=0)
RC4_ARC = RC4_ARC['botflx'].isel(time=0)

RC4_GR = RC4_GR+100000
RC4_ARC = RC4_ARC+100000

lim = 5000000

vmin = 0
vmax = 1

cmap='Blues'
levels = 30

fig, axs = plt.subplots(1,2, figsize=(12, 5), subplot_kw={'projection':ccrs.NorthPolarStereo()},dpi=800)

axs[0].coastlines(resolution='110m',alpha=0.5)
axs[1].coastlines(resolution='110m',alpha=0.5)

RC4_GR.plot(ax=axs[0], transform=ccrs.epsg(6931),cmap=cmap,vmin=vmin,vmax=vmax,levels=levels,add_colorbar=False,alpha=0.5)
RC4_ARC.plot(ax=axs[1], transform=ccrs.epsg(6931),cmap=cmap,vmin=vmin,vmax=vmax,levels=levels,add_colorbar=False,alpha=0.5)


axs[0].set_title('Greenland domain')
axs[1].set_title('Arctic domain')

axs[0].set_ylim([-lim,lim])
axs[0].set_xlim([-lim,lim])

axs[1].set_ylim([-lim,lim])
axs[1].set_xlim([-lim,lim])

plt.tight_layout()

plt.savefig(fig_save_directory+'RACMO_DOMAIN_MAP.png',dpi=800)
