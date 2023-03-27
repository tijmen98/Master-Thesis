#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 13:59:46 2023

@author: tijmen
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
from PIL import Image
from ease_lonlat import EASE2GRID, SUPPORTED_GRIDS

os.chdir('/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Python_scripts')

import Thesis_Functions.calculations as calculations
import Thesis_Functions.data as data
import Thesis_Functions.plotting as plotting

remapdir = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Data/Remap/'

#%%
"""apply MODIS mask to racmo"""



snow_doy_racmo = xr.open_dataset(remapdir+'2001_racmo_snow_doy.nc')
snow_doy_measure = xr.open_dataset(remapdir+'2001_measure_snow_doy.nc')



#%%
fig, axs = plt.subplots(2,2)

axs[0,0].contourf(snow_doy_racmo_masked['first_day'].values)
axs[0,0].set_xlim(250,500)
axs[0,0].set_ylim(250,500)

axs[0,1].contourf(snow_doy_racmo_masked['last_day'].values)
axs[0,1].set_xlim(250,500)
axs[0,1].set_ylim(250,500)

axs[1,0].contourf(snow_doy_racmo['first_day'].values)
axs[1,0].set_xlim(250,500)
axs[1,0].set_ylim(250,500)

axs[1,1].contourf(snow_doy_racmo['last_day'].values)
axs[1,1].set_xlim(250,500)
axs[1,1].set_ylim(250,500)

#%%

fig, axs = plt.subplots(2,2)

axs[0,0].contourf(snow_doy_measure['first_day'].values)
axs[0,0].set_xlim(250,500)
axs[0,0].set_ylim(250,500)

axs[0,1].contourf(snow_doy_measure['last_day'].values)
axs[0,1].set_xlim(250,500)
axs[0,1].set_ylim(250,500)

axs[1,0].contourf(snow_doy_racmo['first_day'].values)
axs[1,0].set_xlim(250,500)
axs[1,0].set_ylim(250,500)

axs[1,1].contourf(snow_doy_racmo['last_day'].values)
axs[1,1].set_xlim(250,500)
axs[1,1].set_ylim(250,500)

plt.show()
#%%
plt.scatter(snow_doy_measure['first_day'][250:500,250:500].values,snow_doy_racmo['first_day'][250:500,250:500].values,s=0.1,zorder=5,color='orange')
plt.xlabel('Satellite')
plt.ylabel('Racmo')
plt.plot(range(365),range(365),color='black',linestyle=(0,(3,3)),zorder=10,alpha=0.5)
plt.xlim(0,365)
plt.ylim(0,365)
ax = plt.gca()
ax.set_aspect(1)
plt.xticks(np.arange(0,365,50))
plt.yticks(np.arange(0,365,50))
plt.tight_layout()