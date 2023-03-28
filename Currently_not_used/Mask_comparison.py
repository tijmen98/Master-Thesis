#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 11:39:28 2023

@author: tijmen
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 08:53:26 2023

@author: tijmen
"""

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from PIL import Image
from ease_lonlat import EASE2GRID, SUPPORTED_GRIDS
import matplotlib as mpl
import cartopy.crs as ccrs

os.chdir('/')

import Thesis_Functions.calculations as calculations
import Thesis_Functions.data as data
import Thesis_Functions.plotting as plotting

remapdir = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Data/Remap/'
fig_save_directory = '/Users/tijmen/Desktop/Figures_Thesis/'

land_ice_mask = xr.DataArray(Image.open(remapdir+'EASE2_N25km.LOCImask_land50_coast0km.720x720.png'))
racmo_ease = xr.open_dataset(remapdir+'RACMO2.4_Snowextend_2001-01-01_EASE_grid.nc')
measure_ease = xr.open_dataset(remapdir+'Measure_2001_EASE.nc')
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


"""Universel land mask"""

land_mask_universal = xr.DataArray(land_mask.to_numpy() * land_mask_racmo.squeeze().isel(time=300).to_numpy())
land_mask_universal = land_mask_universal.rename({'dim_0':'rows','dim_1':'cols'})