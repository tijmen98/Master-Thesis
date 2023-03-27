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

os.chdir('/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Python_scripts')

import Thesis_Functions.calculations as calculations
import Thesis_Functions.data as data
import Thesis_Functions.plotting as plotting

#%%

directory4_DEFAULT = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Testdata/Racmo_2.4/Daily/NC_DEFAULT/'
directory4_MD = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Testdata/Racmo_2.4/Daily/NC_MD/'
remapdir = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Data/Remap/'


RC4_SENSIBLE = (data.Variable_Import(directory4_MD,'hfss'))
RC4_LATENT = (data.Variable_Import(directory4_MD,'hfls'))
RC4_LW_D = (data.Variable_Import(directory4_MD,'rlds'))
RC4_LW_U = -(data.Variable_Import(directory4_MD,'rlus'))
RC4_SW_D = (data.Variable_Import(directory4_MD,'rsds'))
RC4_SW_U = -(data.Variable_Import(directory4_MD,'rsus'))
RC4_GRND = (data.Variable_Import(directory4_MD,'botflx'))
RC4_RFR = (data.Variable_Import(directory4_MD,'rfrzgl'))
RC4_SNM = data.Variable_Import(directory4_MD,'snm')
RC4_SNWD = data.Variable_Import(directory4_DEFAULT,'snd')
RC4_ALB = data.Variable_Import(directory4_DEFAULT,'sswalb')
RC4_TILE = data.Variable_Import(directory4_DEFAULT,'tilefrac8')
RC4_TAS = data.Variable_Import(directory4_DEFAULT,'tas')

t_start = RC4_GRND['time'].values[0]
t_stop = RC4_GRND['time'].values[-1]

directory3_DEFAULT = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Testdata/Racmo_2.3/Daily/NC_DEFAULT/'
directory3_MD = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Testdata/Racmo_2.3/Daily/NC_MD/'

RC3_SENSIBLE = (data.Variable_Import(directory3_MD,'senf')).sel(time=slice(t_start,t_stop))
RC3_LATENT = (data.Variable_Import(directory3_MD,'latf')).sel(time=slice(t_start,t_stop))
RC3_LW_D = (data.Variable_Import(directory3_MD,'lwsd')).sel(time=slice(t_start,t_stop))
RC3_LW_N = (data.Variable_Import(directory3_MD,'lwsn')).sel(time=slice(t_start,t_stop))
RC3_LW_U = RC3_LW_N - RC3_LW_D
RC3_SNM = data.Variable_Import(directory3_MD,'snowmelt').sel(time=slice(t_start,t_stop))


RC3_SW_D = (data.Variable_Import(directory3_MD,'swsd')).sel(time=slice(t_start,t_stop))
RC3_SW_U = (data.Variable_Import(directory3_MD,'swsu')).sel(time=slice(t_start,t_stop))
#RC3_GRND = (data.Variable_Import(directory3_MD,'')).sel(time=slice(t_start,t_stop))
#RC3_RFR = (data.Variable_Import(directory3_MD,'rfrzgl')).sel(time=slice(t_start,t_stop))
RC3_SNWM = data.Variable_Import(directory3_DEFAULT,'snowmass').sel(time=slice(t_start,t_stop))
RC3_SNWRHO = data.Variable_Import(directory3_DEFAULT,'snowrho').sel(time=slice(t_start,t_stop))
RC3_SNWD = (RC3_SNWM*1000)/RC3_SNWRHO
RC3_ALB = data.Variable_Import(directory3_DEFAULT,'alb').sel(time=slice(t_start,t_stop))
RC3_TILE = data.Variable_Import(directory3_DEFAULT,'tileFR8').sel(time=slice(t_start,t_stop))
RC3_TAS = data.Variable_Import(directory3_DEFAULT,'t2m').sel(time=slice(t_start,t_stop))

t_start = RC4_GRND['time'].values[0]
t_stop = RC4_GRND['time'].values[-1]

#RC4_TILE = data.Variable_Import(directory4_DEFAULT,'tilefrac3')
#RC3_TILE = data.Variable_Import(directory3_DEFAULT,'tileFR3').sel(time=slice(t_start,t_stop))

MEASURE_snow_cover = xr.open_dataset(remapdir+'measure_2001_MA.nc')
RACMO_snow_cover = xr.open_dataset(remapdir+'racmo_2001_MA.nc')


RC3_SNWD = (RC3_SNWM*1000)/RC3_SNWRHO

RC4_SNM = data.Variable_Import(directory4_MD,'snm')
RC3_SNM = data.Variable_Import(directory3_MD,'snowmelt').sel(time=slice(t_start,t_stop))


#%%

save_directory = '/Users/tijmen/Desktop/'

areas_int = [['Canada',66.11,-72.84,78],
             ['Greenland-West',68.32,-52.1,59],
             ['Greenland-East',71.15,-23.47,59],
             ['Iceland',65.87,-22.01,60],
             ['Svalbard',77.99,22.77,83]
             ]

i = 4

index = calculations.return_index_from_coordinates(areas_int[i][1],areas_int[i][2],data.Variable_Import(directory4_MD, 'rlds'))

t_start = '2001-03-01'
t_stop = '2001-08-30'

rlat = index[1]
rlon = index[0]


fig, axs = plt.subplots(2,3, figsize=(40,20))
fig.suptitle(areas_int[i][0],fontsize=20)

"""RACMO 4 PLOTTING"""

axs[0,0].plot(RC4_LW_D.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='orange',label='Down')
axs[0,0].plot(RC4_LW_U.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='blue',label='Up')
axs[0,0].plot((RC4_LW_U+RC4_LW_D).isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='black',label='net')
axs[0,0].set_title('Longwave')



axs[1,0].plot(RC4_SW_D.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='orange',label='Down')
axs[1,0].plot(RC4_SW_U.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='blue',label='Up')
axs[1,0].plot((RC4_SW_U+RC4_SW_D).isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='black',label='net')
axs[1,0].set_title('Shortwave')


axs[0,1].plot(RC4_SNWD.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='Blue',label='Snowheight')
axs2=axs[0,1].twinx()
axs2.plot(RC4_SNM.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='Orange',label='Snowmelt')
axs[0,1].set_title('Snowheight')


axs[1,1].plot(RC4_ALB.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='maroon',label='Albedo')
axs[1,1].set_title('Albedo')

axs[1,2].plot(RC4_TAS.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='grey',label='temperature')
axs[1,2].set_title('Temperature')

"""RACMO 3 PLOTTING"""

axs[0,0].plot(RC3_LW_D.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='orange',label='Down',linestyle=(0,(3,3)))
axs[0,0].plot(RC3_LW_U.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='blue',label='Up',linestyle=(0,(3,3)))
axs[0,0].plot((RC3_LW_U+RC3_LW_D).isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='black',label='net',linestyle=(0,(3,3)))
axs[0,0].set_title('Longwave')
axs[0,0].legend()
xlims = axs[0,0].get_xlim()
axs[0,0].hlines(0,xlims[0],xlims[1],color='black',linestyle=(0,(3,3)))

axs[1,0].plot(RC3_SW_D.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='orange',label='Down',linestyle=(0,(3,3)))
axs[1,0].plot(RC3_SW_U.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='blue',label='Up',linestyle=(0,(3,3)))
axs[1,0].plot((RC3_SW_U+RC3_SW_D).isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='black',label='net',linestyle=(0,(3,3)))
axs[1,0].set_title('Shortwave')
axs[1,0].legend()
xlims = axs[0,0].get_xlim()
axs[1,0].hlines(0,xlims[0],xlims[1],color='black',linestyle=(0,(3,3)))

axs[0,1].plot(RC3_SNWD.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='Blue',label='Snowheight',linestyle=(0,(3,3)))
axs2.plot(RC3_SNM.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='Orange',label='Snowmelt',linestyle=(0,(3,3)))
axs2.set_ylabel('kg m-2 s-1')
axs[0,1].set_title('Snow')
axs[0,1].legend()

axs[1,1].plot(RC3_ALB.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='maroon',label='Albedo',linestyle=(0,(3,3)))
axs[1,1].set_title('Albedo')
axs[1,1].legend()

axs[1,2].plot(RC3_TAS.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='grey',label='temperature',linestyle=(0,(3,3)))
axs[1,2].set_title('Temperature')
axs[1,2].legend()
xlims = axs[1,2].get_xlim()
axs[1,2].hlines(273.15,xlims[0],xlims[1],color='black',linestyle=(0,(3,3)))

"""EASE GRID PLOTTING"""

col , row = calculations.return_index_from_coordinates(areas_int[i][1], areas_int[i][2], RACMO_snow_cover)

axs[0,2].plot(MEASURE_snow_cover['merged_snow_cover_extent_moving_average'].isel(cols=col,rows=row).sel(time=slice(t_start,t_stop)),color='Green',label='Satellite Data')
axs[0,2].plot(RACMO_snow_cover['tilefrac5_moving_average'].isel(cols=col,rows=row).sel(time=slice(t_start,t_stop)),color='Orange',label='Racmo')
axs[0,2].legend()


"""apply to all pannels"""

for y in range(len(axs[0,:])):
    for x in range(len(axs[:,0])):
        ylims = axs[x,y].get_ylim()
        axs[x,y].vlines(areas_int[i][3],ylims[0],ylims[1],color='black',linestyle=(0,(3,3)))
        
plt.savefig(save_directory+areas_int[i][0]+'_energyflux_plot')

#%%
   

ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines(alpha=0.5,resolution='50m')


ax.scatter(areas_int[0][2],areas_int[0][1],zorder=10,s=30,label=areas_int[0][0])
ax.scatter(areas_int[1][2],areas_int[1][1],zorder=10,s=30,label=areas_int[1][0])
ax.scatter(areas_int[2][2],areas_int[2][1],zorder=10,s=30,label=areas_int[2][0])
ax.scatter(areas_int[3][2],areas_int[3][1],zorder=10,s=30,label=areas_int[3][0])
ax.scatter(areas_int[4][2],areas_int[4][1],zorder=10,s=30,label=areas_int[4][0])

plt.xlim((-180,180))
plt.ylim((0,90))
plt.legend()
plt.show()     
        
#%%
