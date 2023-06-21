#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 14:33:10 2023

@author: root
"""

import os
import numpy as np
import xarray as xr
import math
import matplotlib.pyplot as plt
import netCDF4 as nc
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as colors
import seaborn as sns
from matplotlib.colors import ListedColormap
import scipy.signal
from matplotlib import gridspec
from matplotlib.lines import Line2D
from PIL import Image
from ease_lonlat import EASE2GRID, SUPPORTED_GRIDS


os.chdir('/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Python_scripts')

import Thesis_Functions.calculations as calculations
import Thesis_Functions.data as data
import Thesis_Functions.plotting as plotting

datadirectory = '/Volumes/Tijmen/Master-Thesis/Data/Snow_cover_analyses/Snow_cover_ease_grid/'
remapdir = '/Volumes/Tijmen/Master-Thesis/Data/Remap/'
year = str(2002)

ds_racmo_acc = xr.open_dataset(datadirectory+year+'/racmo_acc_season.nc')
ds_measure_acc = xr.open_dataset(datadirectory+year+'/measure_acc_season.nc')
da_racmo_acc = ds_racmo_acc['Snowextend Racmo'].squeeze()
da_measure_acc = ds_measure_acc['merged_snow_cover_extent']

ds_racmo_melt = xr.open_dataset(datadirectory+year+'/racmo_melt_season.nc')
ds_measure_melt = xr.open_dataset(datadirectory+year+'/measure_melt_season.nc')
da_racmo_melt = ds_racmo_melt['Snowextend Racmo'].squeeze()
da_measure_melt = ds_measure_melt['merged_snow_cover_extent']


"""Create masks"""

land_ice_mask = xr.DataArray(Image.open(remapdir+'EASE2_N25km.LOCImask_land50_coast0km.720x720.png'))
land_mask_racmo = xr.open_dataset(remapdir+'RACMO2.4_Sea_2001-01-01_EASE_grid.nc')

"""create masks"""

land_mask_racmo = land_mask_racmo.to_array()
land_ice_mask = land_ice_mask[:,:,0]+land_ice_mask[:,:,1]+land_ice_mask[:,:,2]
land_mask = land_ice_mask.where(land_ice_mask>200,0)
land_mask = land_mask.where(land_ice_mask<200,1)

ice_mask = land_ice_mask.where(land_ice_mask<200,0)
ice_mask = ice_mask.where(ice_mask>185,0)
ice_mask = ice_mask.where(ice_mask<175,1)
ice_mask = ice_mask.where(ice_mask==1,np.nan)

sea_mask = land_ice_mask.where(land_ice_mask>150,1)
sea_mask = sea_mask.where(sea_mask<150,0)

ice_mask = ice_mask.rename({'dim_0':'rows','dim_1':'cols'})


"""Universel land mask"""

land_mask_universal = xr.DataArray(land_mask.to_numpy() * land_mask_racmo.squeeze().isel(time=300).to_numpy())
land_mask_universal = land_mask_universal.rename({'dim_0':'rows','dim_1':'cols'})

"""plot season lengths"""

cols = measure_ease.coords['cols'].values
rows = measure_ease.coords['rows'].values

snow_season_length_racmo = xr.open_dataset(savedir_threshold + 'racmo_snow_season_length.nc')
snow_season_length_measure = xr.open_dataset(savedir_threshold + 'measure_snow_season_lenght.nc')

snow_season_length_racmo = snow_season_length_racmo.assign_coords(
    {'latitude': (('rows', 'cols'), latitudes), 'longitude': (('rows', 'cols'), longitudes),
     'rows': (('rows'), rows), 'cols': (('cols'), cols)})
snow_season_length_measure = snow_season_length_measure.assign_coords(
    {'latitude': (('rows', 'cols'), latitudes), 'longitude': (('rows', 'cols'), longitudes),
     'rows': (('rows'), rows), 'cols': (('cols'), cols)})

"""Masking"""

snow_season_length_measure['first_season_length_masked'] = snow_season_length_measure['first_season_length'].where(
    land_mask_universal == 1, np.nan)
snow_season_length_measure['last_season_length_masked'] = snow_season_length_measure['last_season_length'].where(
    land_mask_universal == 1, np.nan)

snow_season_length_racmo['first_season_length_masked'] = snow_season_length_racmo['first_season_length'].where(
    land_mask_universal == 1, np.nan)
snow_season_length_racmo['last_season_length_masked'] = snow_season_length_racmo['last_season_length'].where(
    land_mask_universal == 1, np.nan)

lim = 4000000

vmin = 0
vmax = 150

cmap = 'Blues'
levels = 30

"""First Season"""

fig, axs = plt.subplots(1, 2, figsize=(12, 5), subplot_kw={'projection': ccrs.NorthPolarStereo()}, dpi=800)
fig.suptitle('First season')
axs[0].coastlines(resolution='110m', alpha=0.5)
axs[1].coastlines(resolution='110m', alpha=0.5)

snow_season_length_racmo['first_season_length_masked'].plot(ax=axs[0], transform=ccrs.epsg(6931), cmap=cmap,
                                                            vmin=vmin, vmax=vmax, levels=levels, add_colorbar=False)
snow_season_length_measure['first_season_length_masked'].plot(ax=axs[1], transform=ccrs.epsg(6931), cmap=cmap,
                                                              vmin=vmin, vmax=vmax, levels=levels)

ice_mask.plot(ax=axs[0], cmap='Oranges', add_colorbar=False)
ice_mask.plot(ax=axs[1], cmap='Oranges', add_colorbar=False)

axs[0].set_title('Racmo')
axs[1].set_title('Satellite')

axs[0].set_ylim([-lim, lim])
axs[0].set_xlim([-lim, lim])

axs[1].set_ylim([-lim, lim])
axs[1].set_xlim([-lim, lim])

plt.tight_layout()

plt.savefig(fig_save_directory + 'RACMO_MEASURE_FIRST_SEASON_LENGTH_MAP.png', dpi=800)

"""Last Season"""

vmin = 0
vmax = 100

cmap = 'Blues'
levels = 30

fig, axs = plt.subplots(1, 2, figsize=(12, 5), subplot_kw={'projection': ccrs.NorthPolarStereo()}, dpi=800)
fig.suptitle('Last season')

axs[0].coastlines(resolution='110m', alpha=0.5)
axs[1].coastlines(resolution='110m', alpha=0.5)

snow_season_length_racmo['last_season_length_masked'].plot(ax=axs[0], transform=ccrs.epsg(6931), cmap=cmap,
                                                           vmin=vmin, vmax=vmax, levels=levels, add_colorbar=False)
snow_season_length_measure['last_season_length_masked'].plot(ax=axs[1], transform=ccrs.epsg(6931), cmap=cmap,
                                                             vmin=vmin, vmax=vmax, levels=levels)

ice_mask.plot(ax=axs[0], cmap='Oranges', add_colorbar=False)
ice_mask.plot(ax=axs[1], cmap='Oranges', add_colorbar=False)

axs[0].set_title('Racmo')
axs[1].set_title('Satellite')

axs[0].set_ylim([-lim, lim])
axs[0].set_xlim([-lim, lim])

axs[1].set_ylim([-lim, lim])
axs[1].set_xlim([-lim, lim])

plt.tight_layout()

plt.savefig(fig_save_directory + 'RACMO_MEASURE_LAST_SEASON_LENGTH_MAP.png', dpi=800)

"""Season difference"""

vmin = -50
vmax = 50

cmap = 'seismic'
levels = 30

fig, axs = plt.subplots(1, 2, figsize=(12, 5), subplot_kw={'projection': ccrs.NorthPolarStereo()}, dpi=800)

axs[0].coastlines(resolution='110m', alpha=0.5)
axs[1].coastlines(resolution='110m', alpha=0.5)

(snow_season_length_racmo['first_season_length_masked'] - snow_season_length_measure[
    'first_season_length_masked']).plot(ax=axs[0], transform=ccrs.epsg(6931), cmap=cmap, vmin=vmin, vmax=vmax,
                                        levels=levels, add_colorbar=False)
(snow_season_length_racmo['last_season_length_masked'] - snow_season_length_measure[
    'last_season_length_masked']).plot(ax=axs[1], transform=ccrs.epsg(6931), cmap=cmap, vmin=vmin, vmax=vmax,
                                       levels=levels)

ice_mask.plot(ax=axs[0], cmap='Oranges', add_colorbar=False)
ice_mask.plot(ax=axs[1], cmap='Oranges', add_colorbar=False)

axs[0].set_title('First season')
axs[1].set_title('Last Season')

axs[0].set_ylim([-lim, lim])
axs[0].set_xlim([-lim, lim])

axs[1].set_ylim([-lim, lim])
axs[1].set_xlim([-lim, lim])

plt.tight_layout()

plt.savefig(fig_save_directory + 'RACMO_MEASURE_SEASON_LENGTH_DIFFERENCE_MAP.png', dpi=800)

"""plot scatter heatmap day of year"""