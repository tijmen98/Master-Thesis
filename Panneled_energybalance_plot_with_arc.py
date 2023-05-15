#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 11:39:56 2023

@author: tijmen
"""

import os
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

laptop = True
desktop = False

os.chdir('/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Python_scripts')

import Thesis_Functions.calculations as calculations
import Thesis_Functions.data as data
import Thesis_Functions.plotting as plotting

#%%

if laptop:
    print('Directory structure: laptop')
    datadir = "/Volumes/Tijmen/Master-Thesis/Data/"
if desktop:
    print('Directory structure: desktop')
    datadir = "E:/Master-Thesis/Data/"

racmo_arctic_data_directory = datadir+'RACMO_2.4/PXARC11/2001/'

directory4_DEFAULT = racmo_arctic_data_directory + 'NC_DEFAULT/'
directory4_MD = racmo_arctic_data_directory + 'NC_MD/'
remapdir = datadir+'Remap/'


data.Create_Directory_Information(directory4_MD,'.')
data.Create_Directory_Information(directory4_DEFAULT,'.')

RC4_SENSIBLE = (data.Variable_Import(directory4_MD,'hfss'))
RC4_LATENT = (data.Variable_Import(directory4_MD,'hfls'))
RC4_LW_D = (data.Variable_Import(directory4_MD,'rlds'))
RC4_LW_U = -(data.Variable_Import(directory4_MD,'rlus'))
RC4_SW_D = (data.Variable_Import(directory4_MD,'rsds'))
RC4_SW_U = -(data.Variable_Import(directory4_MD,'rsus'))
RC4_GRND = (data.Variable_Import(directory4_MD,'botflx'))
RC4_RFR = (data.Variable_Import(directory4_MD,'rfrzgl'))
RC4_SNWD = data.Variable_Import(directory4_DEFAULT,'sndp')
RC4_ALB = data.Variable_Import(directory4_DEFAULT,'albcsb')
RC4_TILE = data.Variable_Import(directory4_DEFAULT,'tilefrac8')
RC4_TAS = data.Variable_Import(directory4_DEFAULT,'tas')

year = str(2001)

t_start = year+'-01-01'
t_stop = year+'-12-31'

#%%

save_directory = '/Users/tijmen/Desktop/Figures_Thesis/'

areas_int = [['Canada',66.11,-72.84,78],
             ['Greenland-West',68.32,-52.1,59],
             ['Greenland-East',71.15,-23.47,59],
             ['Iceland',65.87,-22.01,60],
             ['Svalbard',77.99,22.77,83]
             ]

i = 1

index = calculations.return_index_from_coordinates(areas_int[i][1],areas_int[i][2],data.Variable_Import(directory4_MD, 'rlds'))

t_start = '2001-03-01'
t_stop = '2001-08-30'

rlat = index[1]
rlon = index[0]


fig, axs = plt.subplots(2, 3, figsize=(20, 10))
fig.suptitle(areas_int[i][0], fontsize=20)

"""RACMO 4 PLOTTING"""

axs[0,0].plot(RC4_LW_D.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='orange',label='Down')
axs[0,0].plot(RC4_LW_U.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='blue',label='Up')
axs[0,0].plot((RC4_LW_U+RC4_LW_D).isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='black',label='net')
axs[0,0].set_title('Longwave')



axs[1,0].plot(RC4_SW_D.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='orange',label='Down')
axs[1,0].plot(RC4_SW_U.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='blue',label='Up')
axs[1,0].plot((RC4_SW_U+RC4_SW_D).isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='black',label='net')
axs[1,0].set_title('Shortwave')


axs[0,1].plot(RC4_SNWD.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='Blue',label='Racmo')
axs2=axs[0,1].twinx()
axs[0,1].set_title('Snowheight')


axs[1,1].plot(RC4_ALB.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='maroon',label='Racmo')
axs[1,1].set_title('Albedo')

axs[1,2].plot(RC4_TAS.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='grey',label='Racmo')
axs[1,2].set_title('Temperature')


"""apply to all pannels"""

for y in range(len(axs[0,:])):
    for x in range(len(axs[:,0])):
        ylims = axs[x,y].get_ylim()
        axs[x,y].vlines(areas_int[i][3],ylims[0],ylims[1],color='black',linestyle=(0,(3,3)))
        
plt.savefig(save_directory+areas_int[i][0]+'_energyflux_plot')

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
