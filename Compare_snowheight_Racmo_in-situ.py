#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 09:46:58 2023

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
from pyproj import Transformer
from pyproj import CRS
import pandas as pd
import scipy.stats

os.chdir('/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Python_scripts')

import Thesis_Functions.calculations as Calculations
import Thesis_Functions.data as Data
import Thesis_Functions.plotting as Plotting

norway=True
alaska = False

year = 2001



crs_racmo = CRS.from_proj4("-m 57.295779506 +proj=ob_tran +o_proj=latlon +o_lat_p=6.6 +lon_0=180.0")
crs_stations = CRS.from_string("EPSG:4326")

year = str(year)

datadir_in_situ = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Data/In_situ_data/'+year+'/Calculated/'
datadir_racmo = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Data/RACMO_2.4/PXARC11/'

racmo_melt_season_acc = xr.open_dataset(datadir_racmo+'melt_season_snowheight_accumulated.nc')
racmo_acc_season_acc = xr.open_dataset(datadir_racmo+'accumulation_season_snowheight_accumulated.nc')

station_stats = pd.read_csv(datadir_in_situ+'station_statistics_'+year+'.csv',index_col=0)
stationarray = pd.read_csv(datadir_in_situ+'stations_dayly_snowheight_interpolated_'+year+'.csv',index_col=0)

stations = station_stats.columns

fig,axs = plt.subplots(1,1,figsize=(12,12))

locdata = pd.DataFrame(index = ['x','y','in_situ','racmo'],columns=stations)

norm = mpl.colors.LogNorm(vmin=10, vmax=1200)
cmap = mpl.colormaps['Blues']

racmo_melt_season_acc['sndp'].plot(ax=axs,zorder=10,norm=norm, cmap=cmap)



for i,v in enumerate(stations):
    lat = float(station_stats.loc['latitude',v])
    lon = float(station_stats.loc['longitude',v])
    
    y,x = Calculations.return_index_from_coordinates(lat,lon,racmo_acc_season_acc)
    
    
    if x != 0 and x!= len(racmo_acc_season_acc['rlat'])-1 and y != 0 and y!= len(racmo_acc_season_acc['rlon'])-1:
        locdata.loc['x',v]=x
        locdata.loc['y',v]=y
        locdata.loc['in_situ',v]= float(station_stats.loc['melt_season_acc_snowheight',v])
        locdata.loc['racmo',v]= float(racmo_acc_season_acc['sndp'].sel(rlat=racmo_acc_season_acc['rlat'][x],rlon=racmo_acc_season_acc['rlon'][y]).values[0])        
        
        plt.scatter(racmo_acc_season_acc['rlon'][y],racmo_acc_season_acc['rlat'][x],zorder=20, color = cmap(norm(float(station_stats.loc['melt_season_acc_snowheight',v])/100)),edgecolor='black')
    
    else:
        locdata.drop(columns=[v],inplace=True)
    
    
    # if i ==10:
    #     break
    
    
if norway == True:
    
    axs.set_ylim(37,22)
    axs.set_xlim(0,-15)

if alaska == True:
    
    axs.set_ylim(-30,-10)
    axs.set_xlim(0,20)
    

#%%

figlim=100

fig, axs = plt.subplots(1,1,figsize=(10,4),dpi=800)

        
in_situ = locdata.loc['in_situ']/100
racmo = locdata.loc['racmo']

hist, xedges, yedges = np.histogram2d(in_situ,racmo,bins=1)
xidx = np.clip(np.digitize(in_situ, xedges), 0, hist.shape[0]-1)
yidx = np.clip(np.digitize(racmo, yedges), 0, hist.shape[1]-1)
c = hist[xidx, yidx]


RMSE = np.sqrt(np.mean(((racmo-in_situ)**2)))

racmo=racmo.to_list()
in_situ=in_situ.to_list()
regres = scipy.stats.linregress(racmo,in_situ)


axs.set_title('Melt season accumulated snow-height 2001')
axs.scatter(in_situ, racmo, c=c, s = 1, cmap=plt.cm.RdYlBu_r,norm=mpl.colors.LogNorm())
axs.set_xlabel('In_situ [accumulated]')
axs.set_ylabel('Racmo [accumulated]')
axs.plot(range(figlim),range(figlim),color='black',linestyle=(0,(3,3)),zorder=10,alpha=0.5)
axs.set_xlim(0,figlim)
axs.set_ylim(0,figlim)
axs.set_aspect(1)
axs.set_xticks(np.linspace(0,figlim,11))
axs.set_yticks(np.linspace(0,figlim,11))
axs.annotate(('RMSE:'+str(np.round(RMSE,1))),xy=(7*figlim/10,figlim/10+figlim/20))
axs.annotate(('Slope:'+str(np.round(regres.slope,3))),xy=(7*figlim/10,figlim/10))
axs.annotate(('CC:'+str(np.round(regres.rvalue,3))),xy=(7*figlim/10,figlim/10-figlim/20))







