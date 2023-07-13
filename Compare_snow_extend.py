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
version = 'v2'
fig_save_directory = '/Users/tijmen/Desktop/RACMO_24_'+version+'_figures'
datadirectory = '/Volumes/Tijmen/Master-Thesis/Data/Snow_cover_analyses/Snow_cover_ease_grid/'
remapdir = '/Volumes/Tijmen/Master-Thesis/Data/Remap/'


for year in ['2002']:

    print('year is '+year)

    measure_ease = xr.open_dataset(remapdir+'Measure_2001_EASE.nc')

    ds_racmo_acc = xr.open_dataset(datadirectory+year+'/racmo_'+version+'_acc_season.nc')
    ds_measure_acc = xr.open_dataset(datadirectory+year+'/measure_acc_season.nc')
    da_racmo_acc = ds_racmo_acc['Snowextend Racmo'].squeeze()
    da_measure_acc = ds_measure_acc['merged_snow_cover_extent']

    ds_racmo_melt = xr.open_dataset(datadirectory+year+'/racmo_'+version+'_melt_season.nc')
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

    """Latitudes, longitudes"""
    """Create grid with latitudes and longitudes in EASE grid"""

    grid = EASE2GRID(name='EASE2_N25km', **SUPPORTED_GRIDS['EASE2_N25km'])

    latitudes = np.zeros((720, 720))
    longitudes = np.zeros((720, 720))

    for row in range(720):
        for col in range(720):
            try:
                longitudes[row, col], latitudes[row, col] = grid.rc2lonlat(col=col, row=row)

            except:
                True

    """Universel land mask"""

    land_mask_universal = xr.DataArray(land_mask.to_numpy() * land_mask_racmo.squeeze().isel(time=300).to_numpy())
    land_mask_universal = land_mask_universal.rename({'dim_0':'rows','dim_1':'cols'})

    """Assign coordinates"""

    cols = measure_ease.coords['cols'].values
    rows = measure_ease.coords['rows'].values

    ds_measure_acc = ds_measure_acc.assign_coords(
        {'latitude': (('rows', 'cols'), latitudes), 'longitude': (('rows', 'cols'), longitudes),
         'rows': (('rows'), rows), 'cols': (('cols'), cols)})
    ds_racmo_acc = ds_racmo_acc.assign_coords(
        {'latitude': (('rows', 'cols'), latitudes), 'longitude': (('rows', 'cols'), longitudes),
         'rows': (('rows'), rows), 'cols': (('cols'), cols)})

    ds_measure_melt = ds_measure_melt.assign_coords(
        {'latitude': (('rows', 'cols'), latitudes), 'longitude': (('rows', 'cols'), longitudes),
         'rows': (('rows'), rows), 'cols': (('cols'), cols)})
    ds_racmo_melt = ds_racmo_melt.assign_coords(
        {'latitude': (('rows', 'cols'), latitudes), 'longitude': (('rows', 'cols'), longitudes),
         'rows': (('rows'), rows), 'cols': (('cols'), cols)})

    ice_mask = ice_mask.assign_coords({'latitude':(('rows','cols'),latitudes),'longitude':(('rows','cols'),longitudes),'rows':(('rows'),rows),'cols':(('cols'),cols)})


    """Masking"""
    ds_racmo_acc['Snowextent_masked'] = ds_racmo_acc['Snowextend Racmo'].where(
        land_mask_universal == 1, np.nan)
    ds_measure_acc['Snowextent_masked'] = ds_measure_acc['merged_snow_cover_extent'].where(
        land_mask_universal == 1, np.nan)

    ds_racmo_melt['Snowextent_masked'] = ds_racmo_melt['Snowextend Racmo'].where(
        land_mask_universal == 1, np.nan)
    ds_measure_melt['Snowextent_masked'] = ds_measure_melt['merged_snow_cover_extent'].where(
        land_mask_universal == 1, np.nan)

    lim = 4000000

    vmin = 0
    vmax = 100

    cmap = 'Blues'
    levels = 30

    """First Season"""

    fig, axs = plt.subplots(1, 2, figsize=(12, 5), subplot_kw={'projection': ccrs.NorthPolarStereo()}, dpi=400)
    fig.suptitle('Accumulation season')

    ice_mask.plot(ax=axs[0], cmap='Oranges', transform=ccrs.epsg(6931), add_colorbar=False)
    ice_mask.plot(ax=axs[1], cmap='Oranges', transform=ccrs.epsg(6931), add_colorbar=False)

    axs[0].coastlines(resolution='110m', alpha=0.5)
    axs[1].coastlines(resolution='110m', alpha=0.5)

    ds_racmo_acc['Snowextent_masked'].plot(ax=axs[0], transform=ccrs.epsg(6931), cmap=cmap, zorder=10,
                                                                vmin=vmin, vmax=vmax, levels=levels, add_colorbar=False)
    ds_measure_acc['Snowextent_masked'].plot(ax=axs[1], transform=ccrs.epsg(6931), cmap=cmap,
                                                                  vmin=vmin, vmax=vmax, levels=levels, cbar_kwargs={'label':'Snow cover duration [days]'})
    axs[0].set_title('Racmo')
    axs[1].set_title('Satellite')

    axs[0].set_ylim([-lim, lim])
    axs[0].set_xlim([-lim, lim])

    axs[1].set_ylim([-lim, lim])
    axs[1].set_xlim([-lim, lim])

    racmo_mean = np.mean(ds_racmo_acc['Snowextent_masked']).values
    measure_mean = np.mean(ds_measure_acc['Snowextent_masked']).values
    axs[0].annotate('Mean: '+str(np.round(racmo_mean, 2)), (-3500000, -3500000))
    axs[1].annotate('Mean: '+str(np.round(measure_mean, 2)), (-3500000, -3500000))


    plt.tight_layout()

    os.makedirs(fig_save_directory+'/'+year, exist_ok=True)

    plt.savefig(fig_save_directory+'/'+year+'/racmo_'+version+'_measure_acc_season.png', dpi=400)
    plt.close()

    """Last Season"""

    vmin = 0
    vmax = 150

    cmap = 'Blues'
    levels = 30

    fig, axs = plt.subplots(1, 2, figsize=(12, 5), subplot_kw={'projection': ccrs.NorthPolarStereo()}, dpi=400)
    fig.suptitle('Melt season')

    axs[0].coastlines(resolution='110m', alpha=0.5)
    axs[1].coastlines(resolution='110m', alpha=0.5)

    ds_racmo_melt['Snowextent_masked'].plot(ax=axs[0], transform=ccrs.epsg(6931), cmap=cmap,
                                                                vmin=vmin, vmax=vmax, levels=levels, add_colorbar=False)
    ds_measure_melt['Snowextent_masked'].plot(ax=axs[1], transform=ccrs.epsg(6931), cmap=cmap,
                                                                  vmin=vmin, vmax=vmax, levels=levels, cbar_kwargs={'label':'Snow cover duration [days]'})

    ice_mask.plot(ax=axs[0], cmap='Oranges', add_colorbar=False)
    ice_mask.plot(ax=axs[1], cmap='Oranges', add_colorbar=False)

    axs[0].set_title('Racmo')
    axs[1].set_title('Satellite')

    axs[0].set_ylim([-lim, lim])
    axs[0].set_xlim([-lim, lim])

    axs[1].set_ylim([-lim, lim])
    axs[1].set_xlim([-lim, lim])

    racmo_mean = np.mean(ds_racmo_melt['Snowextent_masked']).values
    measure_mean = np.mean(ds_measure_melt['Snowextent_masked']).values
    axs[0].annotate('Mean: '+str(np.round(racmo_mean, 2)), (-3500000, -3500000))
    axs[1].annotate('Mean: '+str(np.round(measure_mean, 2)), (-3500000, -3500000))
    plt.tight_layout()

    plt.savefig(fig_save_directory+'/'+year+'/racmo_'+version+'_measure_melt_season.png', dpi=400)
    plt.close()

    """Season difference"""

    acc_bias = np.nanmean(ds_racmo_acc['Snowextent_masked'] - ds_measure_acc['Snowextent_masked'])
    melt_bias = np.nanmean(ds_racmo_melt['Snowextent_masked'] - ds_measure_melt['Snowextent_masked'])

    acc_rmse = np.sqrt(np.nanmean((ds_racmo_acc['Snowextent_masked'] - ds_measure_acc['Snowextent_masked'])**2))
    melt_rmse= np.sqrt(np.nanmean((ds_racmo_melt['Snowextent_masked'] - ds_measure_melt['Snowextent_masked'])**2))

    vmin = -50
    vmax = 50

    cmap = 'seismic'
    levels = 29

    fig, axs = plt.subplots(1, 2, figsize=(12, 5), subplot_kw={'projection': ccrs.NorthPolarStereo()}, dpi=400)

    axs[0].coastlines(resolution='110m', alpha=0.5)
    axs[1].coastlines(resolution='110m', alpha=0.5)

    (ds_racmo_acc['Snowextent_masked'] - ds_measure_acc['Snowextent_masked']).plot(ax=axs[0], transform=ccrs.epsg(6931), cmap=cmap, vmin=vmin, vmax=vmax,
                                            levels=levels, add_colorbar=False)
    (ds_racmo_melt['Snowextent_masked'] - ds_measure_melt['Snowextent_masked']).plot(ax=axs[1], transform=ccrs.epsg(6931), cmap=cmap, vmin=vmin, vmax=vmax,
                                           levels=levels, cbar_kwargs={'label':'Difference [days]'})

    ice_mask.plot(ax=axs[0], cmap='Oranges', add_colorbar=False)
    ice_mask.plot(ax=axs[1], cmap='Oranges', add_colorbar=False)

    axs[0].set_title('Accumulation season')
    axs[1].set_title('Melt season')

    axs[0].set_ylim([-lim, lim])
    axs[0].set_xlim([-lim, lim])

    axs[1].set_ylim([-lim, lim])
    axs[1].set_xlim([-lim, lim])

    axs[0].annotate('Bias: '+str(np.round(acc_bias, 2)), (-3500000, -3500000))
    axs[1].annotate('Bias: '+str(np.round(melt_bias, 2)), (-3500000, -3500000))

    axs[0].annotate('RMSE: '+str(np.round(acc_rmse, 2)), (-2000000, -3500000))
    axs[1].annotate('RMSE: '+str(np.round(melt_rmse, 2)), (-2000000, -3500000))

    plt.tight_layout()

    plt.savefig(fig_save_directory+'/'+year+'/racmo_'+version+'_measure_season_difference.png', dpi=400)
