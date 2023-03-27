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

os.chdir('/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Python_scripts')

import Thesis_Functions.calculations as calculations
import Thesis_Functions.data as datafunc
import Thesis_Functions.plotting as plotting

directory = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Data/Snow_cover/Measure_2001/'


files = os.listdir(directory)

files = np.sort(files)[2::]

datas = []

for file in files:
    datas.append(xr.open_dataset(directory+file))
    

'merge days to one file'

datamerge = xr.concat(datas,'time')

'change bytes data to snow [1] or no snow [0]'

datamerge['merged_snow_cover_extent'] = datamerge['merged_snow_cover_extent'].where(datamerge['merged_snow_cover_extent']<14,np.nan) * 0 +1
datamerge['merged_snow_cover_extent'] = datamerge['merged_snow_cover_extent'].where(datamerge['merged_snow_cover_extent']==1,0)

datamerge['modis_cloud_gap_filled_snow_cover_extent'] = datamerge['modis_cloud_gap_filled_snow_cover_extent'].where(datamerge['modis_cloud_gap_filled_snow_cover_extent']==10,np.nan) * 0 +1
datamerge['modis_cloud_gap_filled_snow_cover_extent'] = datamerge['modis_cloud_gap_filled_snow_cover_extent'].where(datamerge['modis_cloud_gap_filled_snow_cover_extent']==1,0)


datamerge['merged_snow_cover_extent'].isel(time=0).plot()

filename='Measure_2001_merged.nc'

datamerge.to_netcdf(directory+filename, mode = 'w')

#%%

pyproj.transformer.Tranformer()
