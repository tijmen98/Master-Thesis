#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 10:06:15 2023

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
import pandas as pd

os.chdir('/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Python_scripts')

import Thesis_Functions.calculations as calculations
import Thesis_Functions.data as data
import Thesis_Functions.plotting as plotting

datadir_racmo = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Data/RACMO_2.4/PXARC11/NC_DEFAULT/'
savedir_racmo = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Data/RACMO_2.4/PXARC11/'



year = 2001
breakdate = '/07/01'


year = str(year)


racmo_snowdepth = data.Variable_Import(datadir_racmo, 'sndp').sel(time=slice(year+'/01/01',year+'/12/31'))


melt_season_acc = racmo_snowdepth.sel(time=slice(year+'01/01',year+breakdate)).sum(axis=0)
acc_season_acc = racmo_snowdepth.sel(time=slice(year+breakdate,year+'12/31')).sum(axis=0)

melt_season_acc.to_netcdf(savedir_racmo+'melt_season_snowheight_accumulated.nc')
acc_season_acc.to_netcdf(savedir_racmo+'accumulation_season_snowheight_accumulated.nc')

