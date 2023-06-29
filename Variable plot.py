#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 11:37:20 2023

@author: tijmen
"""

"""t2m plot arctic domain"""

laptop = True
desktop = False

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.colors as colors

os.chdir('/')

import Thesis_Functions.calculations as calculations
import Thesis_Functions.data as data
import Thesis_Functions.plotting as plotting

datadir = "/Volumes/Tijmen/Master-Thesis/Data/"

year = str(2002)
time = 100
safevariable = 'ALBEDO'
variable = 'Clear-sky_albedo'

if laptop:
    print('Directory structure: laptop')
    datadir = "/Volumes/Tijmen/Master-Thesis/Data/"
if desktop:
    print('Directory structure: desktop')
    datadir = "E:/Master-Thesis/Data/"

racmo_ds = xr.open_dataset(datadir+'RACMO_2.4/PXARC11/2001_new/NC_MD/Clearsky_albedo_calculated.nc')
modis_ds = xr.open_dataset(datadir+'RACMO_2.4/PXARC11/2001_v2/NC_MD/Clearsky_albedo_calculated.nc')
racmo_albedo = racmo_ds[variable].sel(time=slice(year+'/01/01',year+'/12/31')).isel(time=time).squeeze()
modis_albedo = modis_ds[variable].sel(time=slice(year+'/01/01',year+'/12/31')).isel(time=time).squeeze()

racmo_albedo = racmo_albedo-modis_albedo

rlon = racmo_ds.rlon.values
rlat = racmo_ds.rlat.values

plt.figure(figsize=(15, 5), dpi=400)

data_crs = ccrs.RotatedPole(pole_longitude=racmo_ds.rotated_pole.grid_north_pole_longitude,
                            pole_latitude=racmo_ds.rotated_pole.grid_north_pole_latitude)
levels=20
norm = colors.Normalize(vmin=-0.2, vmax=0.2)

ax1 = plt.subplot(1, 2, 1, projection=ccrs.Stereographic(central_longitude=0., central_latitude=90.) )
result1 = ax1.contourf(rlon, rlat, racmo_albedo, levels=levels, norm=norm, extend='both', cmap='coolwarm', transform=data_crs)
ax1.coastlines(resolution='50m')
plt.colorbar(result1, orientation='horizontal', label='Albedo', extend='both', fraction=0.046, pad=0.04)
ax1.set_title('Racmo clear sky albedo '+year, size='xx-large')

ax2 = plt.subplot(1, 2, 2, projection=ccrs.Stereographic(central_longitude=0., central_latitude=90.) )
result2 = ax2.contourf(rlon, rlat, modis_albedo, levels=levels, norm=norm, extend='both', cmap='coolwarm', transform=data_crs)
ax2.coastlines(resolution='50m')
plt.colorbar(result2, orientation='horizontal', label='Albedo', extend='both', fraction=0.046, pad=0.04)
ax2.set_title('Modis clear sky albedo '+year, size='xx-large')

plt.savefig('/Users/tijmen/Desktop/Figures_Thesis/Modis_V1-V2.jpeg', dpi=400)

plt.close()

print('')