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

os.chdir('/')

import Thesis_Functions.calculations as calculations
import Thesis_Functions.data as data
import Thesis_Functions.plotting as plotting


datadir = "/Volumes/Tijmen/Master-Thesis/Data/"

in_situ_data_directory = datadir+'In_situ_data/'
racmo_arctic_data_directory = datadir+'RACMO_2.4/PXARC11/2001/'
snow_cover_extend_measure_dir = datadir+'Snow_cover_Measure/'
mask_directory = datadir+'Mask/'

download_measure_dir = 'Download_3-4/'

directory = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Data/RACMO_2.4/PXARC11/NC_DEFAULT/'
remapdir = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Data/Remap/'

t_start = year+'-01-01'
t_stop = year+'-12-30'


tilefrac5 = data.Variable_Import(racmo_arctic_data_directory, 'tilefrac5').sel(time=slice(t_start,t_stop))
tilefrac7 = data.Variable_Import(racmo_arctic_data_directory, 'tilefrac7').sel(time=slice(t_start,t_stop))


tilefrac_57 = (tilefrac5+tilefrac7)
tilefrac_57 = tilefrac_57.rename('Snowextend Racmo')

tilefrac_57.to_netcdf(remapdir+'RACMO2.4_Snowextend_'+t_start+'_RP_grid.nc')


