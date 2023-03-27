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
import cartopy.feature as cfe
import matplotlib.colors as colors
import seaborn as sns
from matplotlib.colors import ListedColormap
import scipy.signal
from matplotlib import gridspec
from matplotlib.lines import Line2D

os.chdir('/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Python_scripts')

import Thesis_Functions.calculations as calculations
import Thesis_Functions.data as data
import Thesis_Functions.plotting as plotting

directory = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Data/Snow_cover/'



ds_racmo = xr.open_dataset(directory+'outfile.nc')
ds_measure = xr.open_dataset(directory+'Measure_2001_merged.nc') 



ds_racmo = ds_racmo.where(ds_racmo>0.4,0)
ds_racmo = ds_racmo.where(ds_racmo==0,1)

date = '2001-03-30'

fig, axs = plt.subplots(1,2, figsize=(40,20))

ds_measure['merged_snow_cover_extent'].sel(time=date).plot(ax=axs[0])
ds_racmo['tilefrac5'].sel(time=date).plot(ax=axs[1])
axs[0].set_xlim((-4000000,2000000))
axs[1].set_xlim((-4000000,2000000))


axs[0].set_ylim((-4000000,2000000))
axs[1].set_ylim((-4000000,2000000))
plt.xlabel('RACMO')
plt.ylabel('MEASURE')

#%%

daterange = ['2001-01-30','2001-03-30']

plt.scatter(ds_measure['merged_snow_cover_extent'].sel(time=daterange).mean(dim='time'),ds_racmo['tilefrac5'].sel(time=daterange).mean(dim='time'))
