#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 11:04:58 2023

@author: tijmen
"""

"""Plot yearly station locations"""


import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from PIL import Image
from ease_lonlat import EASE2GRID, SUPPORTED_GRIDS
import matplotlib as mpl
import cartopy.crs as ccrs
import pandas as pd
from pyproj import Transformer
from pyproj import CRS
import pandas as pd
import scipy.stats
import matplotlib.colors as colors

#change working directory to import thesis functions
os.chdir('/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Python_scripts')

import Thesis_Functions.calculations as Calculations
import Thesis_Functions.data as Data
import Thesis_Functions.plotting as Plotting

#year that is used for the plot
year = str(2001)

#Directory containing the data
datadir = "/Volumes/Tijmen/Master-Thesis/Data/"

#Directory containing the measurement station data
in_situ_data_directory = datadir+'In_situ_data/'

#Directory containing the measurement station data, processed by Data_processing.py and selected for the set year
in_situ_data_directory_year_calculated = in_situ_data_directory+year+'/Calculated/'
#Directory containing racmo model output, used for the projection
racmo_arctic_data_directory = datadir+'RACMO_2.4/PXARC11/2001/'
#Directory where figure is saved
fig_save_directory = '/Users/tijmen/Desktop/Figures_Thesis/'

#read file containing only data from stations that are in the Racmo domain, file created by Data_processing.py
stationselect = pd.read_csv(in_situ_data_directory_year_calculated+'station_in_arctic_domain_'+year+'.csv',index_col=0).to_numpy().flatten('F').tolist()
#read racmo snowhieght data
racmo_24_arc_snowheight = xr.open_dataset(racmo_arctic_data_directory+'NC_DEFAULT/sndp.KNMI-2001.PXARC11.RACMO24_1_complete6_UAR_q_noice_khalo6_era5q.DD.nc').sel(time = slice('2004-01-01','2004-12-31'))
#read file containg station information
station_stats = pd.read_csv(in_situ_data_directory_year_calculated+'station_statistics_'+year+'.csv',index_col=0)

stations = station_stats.columns

levels=20
norm = colors.Normalize(vmin=0,vmax=10)

ds = racmo_24_arc_snowheight

rlon = ds.rlon.values
rlat = ds.rlat.values

plt.figure(figsize=(12,8) )

ax = plt.subplot( projection=ccrs.Stereographic(central_longitude=0., central_latitude=90.))
data_crs = ccrs.RotatedPole(pole_longitude=ds.rotated_pole.grid_north_pole_longitude,
                            pole_latitude=ds.rotated_pole.grid_north_pole_latitude)



# plot data (converting flux to mm/day)
result = ax.contourf(rlon, rlat, racmo_24_arc_snowheight['sndp'].squeeze().isel(time=0), levels=levels, norm=norm,
                     extend='both', cmap='Greys', transform=data_crs)

ax.coastlines(resolution='50m')

plt.colorbar(result, orientation='horizontal', label='Snowdepth', extend='both', fraction=0.046, pad=0.04)

for i,v in enumerate(stations):
    lat = float(station_stats.loc['latitude',v])
    lon = float(station_stats.loc['longitude',v])
    
    y,x = Calculations.return_index_from_coordinates(lat,lon,racmo_24_arc_snowheight)
    
    if v in stationselect:
    
        plt.scatter(racmo_24_arc_snowheight['rlon'][y], racmo_24_arc_snowheight['rlat'][x], zorder=20, color = 'green',
                    transform=data_crs,s=2)
    
plt.savefig(fig_save_directory+'station_locations'+year+'.png')

