#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 11:37:20 2023

@author: tijmen
"""

"""t2m plot arctic domain"""



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


datadir_racmo = '/Volumes/Tijmen/Master-Thesis/Data/RACMO_2.4/PXARC11/NC_DEFAULT/'



year = str(2001)

variable = 'tilefrac9'

ds = xr.open_dataset('/Volumes/Tijmen/Master-Thesis/Data/RACMO_2.4/PXARC11/2001/NC_DEFAULT/'+variable+'.KNMI-2001.PXARC11.RACMO24_1_complete6_UAR_q_noice_khalo6_era5q.DD.nc')

racmo_variable = ds[variable].sel(time=slice(year+'/01/01',year+'/12/31')).isel(time=0).squeeze()


rlon = ds.rlon.values
rlat = ds.rlat.values

plt.figure(figsize=(10,10),dpi=400)

ax = plt.subplot( projection=ccrs.Stereographic(central_longitude=0., central_latitude=90.) )
data_crs = ccrs.RotatedPole(pole_longitude=ds.rotated_pole.grid_north_pole_longitude,
                            pole_latitude=ds.rotated_pole.grid_north_pole_latitude)
levels=20
norm = colors.Normalize(vmin=-25,vmax=10)

# plot data (converting flux to mm/day)
result = ax.contourf(rlon, rlat, racmo_variable, levels = levels , norm=norm, extend='both', cmap='coolwarm', transform=data_crs)

ax.coastlines(resolution='50m')

plt.colorbar(result, orientation='horizontal', label='T2m [\N{DEGREE SIGN}C]', extend='both', fraction=0.046, pad=0.04)

ax.set_title('Yearly mean T2m (2001)', size='xx-large')

plt.savefig('/Users/tijmen/Desktop/Figures_Thesis/'+variable+'_arc.jpeg', dpi=400)

print('')