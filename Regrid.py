#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 10:41:04 2023

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
import pyproj
from ease_lonlat import EASE2GRID, SUPPORTED_GRIDS
import ease_lonlat

os.chdir('/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Python_scripts')

import Thesis_Functions.calculations as calculations
import Thesis_Functions.data as data
import Thesis_Functions.plotting as plotting

directory = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Data/Snow_cover/'

ds_racmo = xr.open_dataset(directory+'RACMO2.4_2001-01-01_2001-12-30.nc')
ds_measure = xr.open_dataset(directory+'Measure_2001_merged.nc')

#%%
# or using parameters taken from NSIDC and kept in SUPPORTED_GRIDS
grid = EASE2GRID(name='EASE2_N25km', **SUPPORTED_GRIDS['EASE2_N25km'])

latitudes=np.zeros((720,720))
longitudes=np.zeros((720,720))

for row in range(720):
    for col in range(720):
        try:
            longitudes[row,col],latitudes[row,col] = grid.rc2lonlat(col=col,row=row)
            
        except:
            True
            


ds_measure['latitude']=(('rows','cols'),latitudes)
ds_measure['longitude']=(('rows','cols'),longitudes)
ds_measure['cols']=range(720)
ds_measure['rows']=range(720)
