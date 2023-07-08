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
startmonth = '03'
endmonth= '04'
safevariable = 'ALBEDO'
variable = 'Clear-sky_albedo'

if laptop:
    print('Directory structure: laptop')
    datadir = "/Volumes/Tijmen/Master-Thesis/Data/"
if desktop:
    print('Directory structure: desktop')
    datadir = "E:/Master-Thesis/Data/"

tilefrac6 = xr.open_dataset(datadir+'RACMO_2.4/PXARC11/2001_new/NC_DEFAULT/tilefrac6.KNMI-2001.PXARC11.RACMO24_1_complete6_UAR_q_noice_khalo6_era5q.DD.nc')['tilefrac6']
tilefrac7 = xr.open_dataset(datadir+'RACMO_2.4/PXARC11/2001_new/NC_DEFAULT/tilefrac7.KNMI-2001.PXARC11.RACMO24_1_complete6_UAR_q_noice_khalo6_era5q.DD.nc')['tilefrac7']

forest = tilefrac6.isel(time=0)+tilefrac7.isel(time=0)

forest = forest.where(forest>0.4).squeeze()


var1_ds = xr.open_dataset(datadir+'RACMO_2.4/PXARC11/2001_new/NC_MD/Clearsky_albedo_calculated.nc')
var2_ds = xr.open_dataset(datadir+'RACMO_2.4/PXARC11/2001_v2/NC_MD/Clearsky_albedo_calculated.nc')
var1_dr = var1_ds[variable].sel(time=slice(year+'/'+startmonth+'/01', year+'/'+endmonth+'/01')).mean(dim='time').squeeze()
var2_dr = var2_ds[variable].sel(time=slice(year+'/'+startmonth+'/01', year+'/'+endmonth+'/01')).mean(dim='time').squeeze()

var_dif_dr = var2_dr-var1_dr

rlon = var1_ds.rlon.values
rlat = var1_ds.rlat.values

plt.figure(figsize=(8, 8), dpi=400)

data_crs = ccrs.RotatedPole(pole_longitude=var1_ds.rotated_pole.grid_north_pole_longitude,
                            pole_latitude=var1_ds.rotated_pole.grid_north_pole_latitude)
levels=20
norm = colors.Normalize(vmin=-0.1, vmax=0.1)

ax1 = plt.subplot(1, 1, 1, projection=ccrs.Stereographic(central_longitude=0., central_latitude=90.) )
result1 = ax1.contourf(rlon, rlat, var_dif_dr, levels=levels, norm=norm, extend='both', cmap='coolwarm', transform=data_crs)
ax1.coastlines(resolution='50m')
plt.colorbar(result1, orientation='horizontal', label='Albedo', extend='both', fraction=0.046, pad=0.04)
ax1.set_title('Racmo clear sky albedo April '+year, size='xx-large')

plt.savefig('/Users/tijmen/Desktop/Figures_Thesis/V1-V2_albedo_.jpeg', dpi=230)

plt.close()

print('')