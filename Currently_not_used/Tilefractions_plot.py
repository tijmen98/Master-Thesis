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

import Thesis_Functions.calculations as calculations
import Thesis_Functions.data as data
import Thesis_Functions.plotting as plotting

directory = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/TestData/Racmo 2.4/Daily/NC_DEFAULT/'

data.Create_Directory_Information(directory, '.')

areas_int = [['Canada',66.11,-72.84],
             ['Greenland',68.32,-52.1],
             ['Iceland',65.87,-22.01],
             ['Svalbard',77.99,22.77]
             ]
tilefrac = data.Variable_Import(directory, 'tilefrac1')

#%%

# create list of colors
color_legend = [['Tilefrac 1','Tilefrac 2','Interception','Low vegetation','Low vegetation, covered in snow','High vegetation','High vegetation, covered in snow','Bare ground','Tilefrac 9, lakes','Tilefrac 10'],
                ['white','white','tomato','lawngreen','turquoise','green','blue','brown','purple','orange']]

for i in range(len(areas_int)):
    
    fig, axs = plt.subplots(1,1,figsize=(20,16))
    
    index = calculations.return_index_from_coordinates(areas_int[i][1],areas_int[i][2],tilefrac)
    
    cumulative_fractions = np.zeros([10,len(tilefrac)])
    
    for tile in range(10):
        cumulative_fractions[tile]=data.Variable_Import(directory, 'tilefrac'+str(tile+1)).isel(rlon=index[0],rlat=index[1]).values[:,0]

    plt.stackplot(tilefrac['time'].values,cumulative_fractions,colors=color_legend[1])
    
    plt.title(areas_int[i][0])
    
    plt.legend(color_legend[0])

    plt.show()
            

