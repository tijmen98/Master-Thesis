#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 09:46:58 2023

@author: tijmen
"""

import imageio
import matplotlib as mpl
import pandas as pd
from pyproj import CRS
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.colors as colors
import geopandas as gpd
import cartopy.feature as cfeature

os.chdir('/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Python_scripts')

import Thesis_Functions.calculations as Calculations
import Thesis_Functions.data as Data
import Thesis_Functions.plotting as Plotting

crs_racmo = CRS.from_proj4("-m 57.295779506 +proj=ob_tran +o_proj=latlon +o_lat_p=6.6 +lon_0=180.0")
crs_stations = CRS.from_string("EPSG:4326")


years = [2002]

"""Directories"""

datadir = "/Volumes/Tijmen/Master-Thesis/Data/"

modis_data_directory = datadir+'MODIS/'
racmo_arctic_data_directory = datadir+'RACMO_2.4/PXARC11/2001_new/'
mask_directory = datadir+'Mask/'
fig_save_directory = '/Users/tijmen/Desktop/Figures_Thesis/'

images = []

tilefrac6 = xr.open_dataset(racmo_arctic_data_directory+'NC_DEFAULT/tilefrac6.KNMI-2001.PXARC11.RACMO24_1_complete6_UAR_q_noice_khalo6_era5q.DD.nc')['tilefrac6']
tilefrac7 = xr.open_dataset(racmo_arctic_data_directory+'NC_DEFAULT/tilefrac7.KNMI-2001.PXARC11.RACMO24_1_complete6_UAR_q_noice_khalo6_era5q.DD.nc')['tilefrac7']
forest = tilefrac6.isel(time=0)+tilefrac7.isel(time=0)
forest = forest.where(forest>0.4).squeeze()

rlon= tilefrac6.rlon.values
rlat= tilefrac6.rlat.values

albedo_racmo = xr.open_dataset(racmo_arctic_data_directory+'NC_MD/Clearsky_albedo_calculated.nc')
bounds = np.linspace(0, 1, 20)
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256, extend='both')

for year in years:
    albedo_modis_year = xr.open_dataset(modis_data_directory + year+'_RCG.nc')
    albedo_racmo_year = albedo_racmo.sel(time=slice(year+'/01/01', year+'/12/31'))['Clear-sky_albedo']

    rlon = albedo_racmo_year.rlon.values
    rlat = albedo_racmo_year.rlat.values

    for day in range(len(albedo_racmo_year['time'])):

        albedo_racmo_day = albedo_racmo_year.isel(time=day)
        albedo_modis_day = albedo_modis_year.isel(time=day)

        cmap = 'seismic'
        fig = plt.figure(figsize=(12, 10))
        ax1 = fig.add_subplot(1, 2, 1, projection=ccrs.Stereographic(central_longitude=0., central_latitude=90.) )
        ax2 = fig.add_subplot(1, 2, 2, projection=ccrs.Stereographic(central_longitude=0., central_latitude=90.))

        ax1.add_feature(cfeature.OCEAN)
        ax1.coastlines(resolution='50m')
        ax1.add_feature(cfeature.LAND)

        ax1.contourf(rlon, rlat, albedo_racmo_day, cmap=cmap, transform = ccrs.RotatedPole(pole_longitude=albedo_racmo.rotated_pole.grid_north_pole_longitude,
                                pole_latitude=albedo_racmo.rotated_pole.grid_north_pole_latitude))

        ax2.contourf(rlon, rlat, albedo_modis_day, cmap=cmap, transform = ccrs.RotatedPole(pole_longitude=albedo_racmo.rotated_pole.grid_north_pole_longitude,
                                pole_latitude=albedo_racmo.rotated_pole.grid_north_pole_latitude))

        ax1.set_xlim(-3939996.3398490055, 4199452.516011998)
        ax1.set_ylim(-4171531.5955094383, 4399068.623679058)

        axpos = ax1.get_position()
        pos_x = axpos.x0 + axpos.width + 0.01
        pos_y = axpos.y0
        cax_width = 0.01
        cax_height = axpos.height
        pos_cax = fig.add_axes([pos_x, pos_y, cax_width, cax_height])
        plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=pos_cax, label='Bias [cm]')

        axpos = ax2.get_position()
        pos_x = axpos.x0 + axpos.width + 0.01
        pos_y = axpos.y0
        cax_width = 0.01
        cax_height = axpos.height
        pos_cax = fig.add_axes([pos_x, pos_y, cax_width, cax_height])
        plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=pos_cax, label='Bias [cm]')

        fig.set_title('Albedo year:'+ year + 'day' +str(day), size='xx-large')

        plt.savefig(fig_save_directory+year+'/Albedomovie/'+str(day)+'_albedo_plot.png')

        plt.close(fig)

        images.append(imageio.v3.imread(fig_save_directory+year+'/Albedomovie/'+str(day)+'_albedo_plot.png'))

imageio.mimsave(fig_save_directory+'movie_albedo.gif', images)
