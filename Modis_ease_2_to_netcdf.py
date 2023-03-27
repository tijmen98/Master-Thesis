#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 20:56:04 2023

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
land_ice_mask = xr.DataArray(Image.open(remapdir+'EASE2_N6.25km.LOCImask_land50_coast0km.2880x2880.png'))

"""create masks"""


land_ice_mask_tot = land_ice_mask[:,:,0]+land_ice_mask[:,:,1]+land_ice_mask[:,:,2]

land_mask = land_ice_mask_tot.where(land_ice_mask_tot>200,0)
land_mask = land_mask.where(land_ice_mask_tot<200,1)
ice_mask = land_ice_mask_tot.where(land_ice_mask_tot<200,0)
ice_mask = ice_mask.where(ice_mask>175,0)
ice_mask = ice_mask.where(ice_mask<175,1)
sea_mask = land_ice_mask_tot.where(land_ice_mask_tot>150,1)
sea_mask = sea_mask.where(sea_mask<150,0)


grid = EASE2GRID(name='EASE2_N625km', **SUPPORTED_GRIDS['EASE2_N625km'])

latitudes=np.zeros((2880,2880))
longitudes=np.zeros((2880,2880))

for row in range(2880):
    for col in range(2880):
        try:
            longitudes[row,col],latitudes[row,col] = grid.rc2lonlat(col=col,row=row)
            
        except:
            True

#%%


land_mask = land_ice_mask
land_mask['ice'] = ice_mask
land_mask['land'] = land_mask
land_mask['sea'] = sea_mask
land_mask.rename({'dim_0':'rows','dim_1':'cols'})
