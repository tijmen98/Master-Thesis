#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 12:01:55 2023

@author: tijmen

Determine the DOY that snow has melted from 
"""

import os
import numpy as np
import xarray as xr
import math
import matplotlib.pyplot as plt
import netCDF4
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfe
import matplotlib.colors as colors
import seaborn as sns
from matplotlib.colors import ListedColormap
import scipy.signal
from matplotlib import gridspec
from matplotlib.lines import Line2D
import pyproj

years = [2003, 2004]


os.chdir('/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Python_scripts')

import Thesis_Functions.calculations as calculations
import Thesis_Functions.data as datafunc
import Thesis_Functions.plotting as plotting

"""Directories"""

datadir = "/Volumes/Tijmen/Master-Thesis/Data/"

snow_cover_extend_measure_dir = datadir+'Snow_cover_Measure/'
download_measure_dir = 'Download_3-4/'

filename='Measure_multiple_merged.nc'

directories = os.listdir(snow_cover_extend_measure_dir+download_measure_dir)

directories = sorted((f for f in directories if not f.startswith(".")), key=str.lower)

data = []

for directory in directories:
    files = os.listdir(snow_cover_extend_measure_dir + download_measure_dir + directory)
    file = sorted(file for file in files if not file.startswith(".") and file.endswith(".nc"))

    data.append(xr.open_dataset(snow_cover_extend_measure_dir + download_measure_dir + directory +'/' + file[0]))


'merge days to one file'

datamerge = xr.concat(data,'time')

'change bytes data to snow [1] or no snow [0]'

datamerge['merged_snow_cover_extent'] = datamerge['merged_snow_cover_extent'].where(datamerge['merged_snow_cover_extent']<14,np.nan) * 0 +1
datamerge['merged_snow_cover_extent'] = datamerge['merged_snow_cover_extent'].where(datamerge['merged_snow_cover_extent']==1,0)

datamerge['modis_cloud_gap_filled_snow_cover_extent'] = datamerge['modis_cloud_gap_filled_snow_cover_extent'].where(datamerge['modis_cloud_gap_filled_snow_cover_extent']==10,np.nan) * 0 +1
datamerge['modis_cloud_gap_filled_snow_cover_extent'] = datamerge['modis_cloud_gap_filled_snow_cover_extent'].where(datamerge['modis_cloud_gap_filled_snow_cover_extent']==1,0)


datamerge['merged_snow_cover_extent'].isel(time=900).plot()

datamerge.to_netcdf(snow_cover_extend_measure_dir+filename)
#%%
