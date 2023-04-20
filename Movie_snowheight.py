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


years = [2002, 2003, 2004]
months=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
month_names = ['January', 'February', 'March', ' April', 'May', 'June',
               'July', 'August', 'September', 'October', 'November', 'December']

"""Directories"""

datadir = "/Volumes/Tijmen/Master-Thesis/Data/"

in_situ_data_directory = datadir+'In_situ_data/'
racmo_arctic_data_directory = datadir+'RACMO_2.4/PXARC11/2001/'
snow_cover_extend_measure_dir = datadir+'Snow_cover_Measure/'
mask_directory = datadir+'Mask/'
fig_save_directory = '/Users/tijmen/Desktop/Figures_Thesis/'

images = []

tilefrac6 = xr.open_dataset(racmo_arctic_data_directory+'NC_DEFAULT/tilefrac6.KNMI-2001.PXARC11.RACMO24_1_complete6_UAR_q_noice_khalo6_era5q.DD.nc')['tilefrac6']
tilefrac7 = xr.open_dataset(racmo_arctic_data_directory+'NC_DEFAULT/tilefrac7.KNMI-2001.PXARC11.RACMO24_1_complete6_UAR_q_noice_khalo6_era5q.DD.nc')['tilefrac7']

rlon= tilefrac6.rlon.values
rlat= tilefrac6.rlat.values

forest = tilefrac6.isel(time=0)+tilefrac7.isel(time=0)

forest = forest.where(forest>0.4).squeeze()

for year in years:

    year = str(year)
    print('Starting calculations for year: ' + year)
    """Year specific directories"""

    in_situ_data_directory_year = in_situ_data_directory + year + '/'
    in_situ_data_directory_year_calculated = in_situ_data_directory + year + '/Calculated/'

    ds = xr.open_dataset(racmo_arctic_data_directory + 'NC_DEFAULT/sndp.KNMI-2001.PXARC11.RACMO24_1_complete6_UAR_q_noice_khalo6_era5q.DD.nc')
    racmo_sndp = ds['sndp'].sel(time=slice(year + '/01/01', year + '/12/31')).mean(dim='time').squeeze()
    rlon = ds.rlon.values
    rlat = ds.rlat.values

    station_locs = pd.read_csv(in_situ_data_directory_year_calculated + 'station_in_arctic_domain_' + year + '.csv',
                               index_col=0)

    for month in months:
        print('Month ' + str(month))

        monthdir_in_situ = in_situ_data_directory_year_calculated+'snow_depth/'+ 'month_' + str(month)
        monthdir_racmo = racmo_arctic_data_directory + year + 'month_' + str(month)

        """Import statistics"""

        vmin= -20
        vmax = 20

        station_statistics = pd.read_csv(monthdir_in_situ+'/Calculated_statistics.csv', index_col=0)
        bias = station_statistics.loc['bias']
        racmo_mean = np.round(np.mean(station_statistics.loc['racmo_mean'].values))
        in_situ_mean = np.round(np.mean(np.isnan(station_statistics.loc['in_situ_mean'].values))*100)

        bounds = np.linspace(vmin, vmax, 10)

        """Plot basemap"""

        norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256, extend='both')
        cmap = 'seismic'
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.Stereographic(central_longitude=0., central_latitude=90.) )

        station_lat = np.array(station_locs.loc['latitude'].values).astype(float)
        station_lon = np.array(station_locs.loc['longitude'].values).astype(float)
        gdf = gpd.GeoDataFrame(bias, geometry=gpd.points_from_xy(station_lon, station_lat,
                                                                 crs=ccrs.Stereographic(central_latitude=90.,
                                                                                        central_longitude=0.)))

        ax.add_feature(cfeature.OCEAN)
        ax.coastlines(resolution='50m')
        ax.add_feature(cfeature.LAND)

        ax.contourf(rlon, rlat, forest, cmap='Greens', transform = ccrs.RotatedPole(pole_longitude=ds.rotated_pole.grid_north_pole_longitude,
                            pole_latitude=ds.rotated_pole.grid_north_pole_latitude))

        stat_dots = gdf.plot(column='bias', ax=ax, markersize=20, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())

        ax.set_xlim(-3939996.3398490055, 4199452.516011998)
        ax.set_ylim(-4171531.5955094383, 4399068.623679058)

        axpos = ax.get_position()
        pos_x = axpos.x0 + axpos.width + 0.01
        pos_y = axpos.y0
        cax_width = 0.01
        cax_height = axpos.height
        pos_cax = fig.add_axes([pos_x, pos_y, cax_width, cax_height])
        plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=pos_cax, label='Bias [cm]')

        ax.set_title(month_names[month-1] + ' ' + year, size='xx-large')

        plt.savefig(fig_save_directory+year+'/month_'+str(month)+'_bias_plot.png')

        plt.close(fig)

        images.append(imageio.v3.imread(fig_save_directory+year+'/month_'+str(month)+'_bias_plot.png'))
        images.append(imageio.v3.imread(fig_save_directory + year + '/month_' + str(month) + '_bias_plot.png'))
        images.append(imageio.v3.imread(fig_save_directory + year + '/month_' + str(month) + '_bias_plot.png'))

        del gdf

imageio.mimsave(fig_save_directory+'movie.gif', images)
