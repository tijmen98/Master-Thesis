#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 11:39:56 2023

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

os.chdir('/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Python_scripts')

import Thesis_Functions.calculations as calculations
import Thesis_Functions.data as data
import Thesis_Functions.plotting as plotting

directory = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Data/RACMO_2.4/PXARC11/NC_DEFAULT/'
savedir = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Data/Remap/'
variable1 = 'tilefrac5'
variable2 = 'tilefrac7'
variable3 = 'tilefrac1'
variable4 = 'tilefrac2'

t_start = '2001-01-01'
t_stop = '2001-12-30'


tilefrac5 = data.Variable_Import(directory, variable1).sel(time=slice(t_start,t_stop))
tilefrac7 = data.Variable_Import(directory, variable2).sel(time=slice(t_start,t_stop))

tilefrac1 = data.Variable_Import(directory, variable3).sel(time=slice(t_start,t_stop))
tilefrac2 = data.Variable_Import(directory, variable4).sel(time=slice(t_start,t_stop))

tilefrac_57 = (tilefrac5+tilefrac7)
tilefrac_57 = tilefrac_57.rename('tile')

tilefrac_57.to_netcdf(savedir+'RACMO2.4_Snowextend_'+t_start+'_RP_grid.nc')


sea = ((tilefrac1+tilefrac2)-1)*-1

sea.to_netcdf(savedir+'RACMO2.4_Sea_'+t_start+'_RP_grid.nc')

