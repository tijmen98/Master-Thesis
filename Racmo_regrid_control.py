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
import scipy.stats

os.chdir('/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Python_scripts')

import Thesis_Functions.calculations as calculations
import Thesis_Functions.data as data
import Thesis_Functions.plotting as plotting

threshold = 0.4 #threshold where difference between snow and no snow is classified
year = 2001


datadir = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Data/Snow_cover_analyses/Snow_cover_ease_grid/'
savedir_threshold = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Data/Snow_cover_analyses/Threshold analyses/'+str(year)+'/Threshold_'+str(threshold)+'/'

remapdir = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Data/Remap/'
fig_save_directory = '/Users/tijmen/Desktop/Figures_Thesis/Threshold_'+str(threshold)+'/'

Day_of_year_calculation = False
Season_length_calculation = True
Day_of_year_plot = False
Season_length_plot = True
Location_map = False



land_ice_mask = xr.DataArray(Image.open(remapdir+'EASE2_N25km.LOCImask_land50_coast0km.720x720.png'))
racmo_ease = xr.open_dataset(datadir+'RACMO2.4_Snowextend_2001-01-01_EASE_grid.nc')
measure_ease = xr.open_dataset(datadir+'Measure_2001_EASE.nc')
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


"""Create grid with latitudes and longitudes in EASE grid"""

grid = EASE2GRID(name='EASE2_N25km', **SUPPORTED_GRIDS['EASE2_N25km'])

latitudes=np.zeros((720,720))
longitudes=np.zeros((720,720))

for row in range(720):
    for col in range(720):
        try:
            longitudes[row,col],latitudes[row,col] = grid.rc2lonlat(col=col,row=row)
            
        except:
            True
racmo_ease['lat']=(('rows','cols'),latitudes)
racmo_ease['lon']=(('rows','cols'),longitudes)

measure_ease['lat']=(('rows','cols'),latitudes)
measure_ease['lon']=(('rows','cols'),longitudes)

"""Make racmo data boolian"""



time = '2001-03-20'

racmo_ease['tile'] = racmo_ease['tile'].where(racmo_ease['tile']>0,np.nan)
racmo_ease['tile'] = racmo_ease['tile'].where(racmo_ease['tile']>threshold,np.nan)*0+1
racmo_ease['tile'] = racmo_ease['tile'].where(racmo_ease['tile']==1,0)


racmo_land_masked = land_mask * racmo_ease['tile'].squeeze().sel(time=time).to_numpy()
racmo_sea_masked = sea_mask * racmo_ease['tile'].squeeze().sel(time=time).to_numpy()
racmo_ice_masked = ice_mask * racmo_ease['tile'].squeeze().sel(time=time).to_numpy()

"""get moving average and save to netcdf with interval 'moving_average' """

moving_average = 7 #number of values that are averaged


racmo_ease['tile_moving_average'] = racmo_ease['tile'].rolling(time=moving_average, center=True).mean()
measure_ease['merged_snow_cover_extent_moving_average'] = measure_ease['merged_snow_cover_extent'].rolling(time=moving_average, center=True).mean()


racmo_ease['lat']=(('rows','cols'),latitudes)
racmo_ease['lon']=(('rows','cols'),longitudes)

measure_ease['lat']=(('rows','cols'),latitudes)
measure_ease['lon']=(('rows','cols'),longitudes)

cols = measure_ease.coords['cols'].values
rows = measure_ease.coords['rows'].values


ice_mask = ice_mask.assign_coords({'latitude':(('rows','cols'),latitudes),'longitude':(('rows','cols'),longitudes),'rows':(('rows'),rows),'cols':(('cols'),cols)})



if True == False:
    measure_ease.to_netcdf(savedir_threshold+'measure_MA.nc')
    racmo_ease.to_netcdf(savedir_threshold+'racmo_MA.nc')


"""Get doy of year that snow is gone"""

if Day_of_year_calculation == True:
    rows = measure_ease['rows'].values
    cols = measure_ease['cols'].values
    
    first_day_measure = np.zeros((len(rows),len(cols)))
    last_day_measure = np.zeros((len(rows),len(cols)))
    
    for col in range(len(cols)):
        for row in range(len(rows)):
            no_snow_measure = np.where(measure_ease['merged_snow_cover_extent_moving_average'].isel(cols=col,rows=row).values==0)
            try:
                first_day_measure[row,col] = no_snow_measure[0][0]
                last_day_measure[row,col] = no_snow_measure[0][-1]
            except:
                first_day_measure[row,col] = np.nan
                last_day_measure[row,col] = np.nan       
                
    first_day_racmo = np.zeros((len(rows),len(cols)))
    last_day_racmo = np.zeros((len(rows),len(cols)))
    
    for col in range(len(cols)):
        for row in range(len(rows)):
            no_snow_racmo = np.where(racmo_ease['tile_moving_average'].isel(cols=col,rows=row).values==0)
            try:
                first_day_racmo[row,col] = no_snow_racmo[0][0]
                last_day_racmo[row,col] = no_snow_racmo[0][-1]
            except:
                first_day_racmo[row,col] = np.nan
                last_day_racmo[row,col] = np.nan   
    
    """merge snow data and save to netcdf"""
                
    first_day_measure = xr.DataArray(first_day_measure,name='first_day')
    last_day_measure = xr.DataArray(last_day_measure,name='last_day')
    
    snow_doy_measure = xr.merge([first_day_measure,last_day_measure])
    snow_doy_measure = snow_doy_measure.rename_dims({'dim_0':'rows','dim_1':'cols'})
    
    first_day_racmo = xr.DataArray(first_day_racmo,name='first_day')
    last_day_racmo = xr.DataArray(last_day_racmo,name='last_day')
    
    snow_doy_racmo = xr.merge([first_day_racmo,last_day_racmo])
    snow_doy_racmo = snow_doy_racmo.rename_dims({'dim_0':'rows','dim_1':'cols'})
    
    snow_doy_racmo.to_netcdf(savedir_threshold+'racmo_snow_doy.nc')
    snow_doy_measure.to_netcdf(savedir_threshold+'measure_snow_doy.nc')
    
"""get number of days with snowcover from 01-01 to 'breakdate' and from 'breakdate' to 31-12"""


breakdate = '07/01'
year = str(year)+'/'

"""Calculate season lengths"""

if Season_length_calculation == True:
    rows = measure_ease['rows'].values
    cols = measure_ease['cols'].values
    
    snow_first_measure = np.zeros((len(rows),len(cols)))
    snow_last_measure = np.zeros((len(rows),len(cols)))
    
    for col in range(len(cols)):
        for row in range(len(rows)):
            try:
                snow_first_measure[row,col] = len(np.where(measure_ease['merged_snow_cover_extent'].isel(cols=col,rows=row).sel(time = slice(year+'01/01',year+breakdate)).values==1)[0])
                snow_last_measure[row,col] = len(np.where(measure_ease['merged_snow_cover_extent'].isel(cols=col,rows=row).sel(time = slice(year+breakdate,year+'12-31')).values==1)[0])
            except:
                snow_first_measure[row,col] = np.nan
                snow_first_measure[row,col] = np.nan       
                
    snow_first_racmo = np.zeros((len(rows),len(cols)))
    snow_last_racmo = np.zeros((len(rows),len(cols)))
    
    for col in range(len(cols)):
        for row in range(len(rows)):
            try:
                snow_first_racmo[row,col] = len(np.where(racmo_ease['tile'].isel(cols=col,rows=row).sel(time = slice(year+'01/01',year+breakdate)).values==1)[0])
                snow_last_racmo[row,col] = len(np.where(racmo_ease['tile'].isel(cols=col,rows=row).sel(time = slice(year+breakdate,year+'12/31')).values==1)[0])
            except:
                snow_first_racmo[row,col] = np.nan
                snow_first_racmo[row,col] = np.nan     
    
    """merge snow data and save to netcdf"""
                
    snow_first_measure = xr.DataArray(snow_first_measure,name='first_season_length')
    snow_last_measure = xr.DataArray(snow_last_measure,name='last_season_length')
    
    snow_season_length_measure = xr.merge([snow_first_measure,snow_last_measure])
    snow_season_length_measure = snow_season_length_measure.rename_dims({'dim_0':'rows','dim_1':'cols'})
    
    snow_first_racmo = xr.DataArray(snow_first_racmo,name='first_season_length')
    snow_last_racmo = xr.DataArray(snow_last_racmo,name='last_season_length')
    
    snow_season_length_racmo = xr.merge([snow_first_racmo,snow_last_racmo])
    snow_season_length_racmo = snow_season_length_racmo.rename_dims({'dim_0':'rows','dim_1':'cols'})
    
    snow_season_length_racmo = snow_season_length_racmo.assign_coords({'latitude':(('rows','cols'),latitudes),'longitude':(('rows','cols'),longitudes),'rows':(('rows'),rows),'cols':(('cols'),cols)})
    snow_season_length_measure = snow_season_length_measure.assign_coords({'latitude':(('rows','cols'),latitudes),'longitude':(('rows','cols'),longitudes),'rows':(('rows'),rows),'cols':(('cols'),cols)})

    """Masking"""

    snow_season_length_measure['first_season_length_masked'] = snow_season_length_measure['first_season_length'].where(land_mask_universal==1,np.nan)
    snow_season_length_measure['last_season_length_masked'] = snow_season_length_measure['last_season_length'].where(land_mask_universal==1,np.nan)

    snow_season_length_racmo['first_season_length_masked'] = snow_season_length_racmo['first_season_length'].where(land_mask_universal==1,np.nan)
    snow_season_length_racmo['last_season_length_masked'] = snow_season_length_racmo['last_season_length'].where(land_mask_universal==1,np.nan)
    
    snow_season_length_racmo.to_netcdf(savedir_threshold+'racmo_snow_season_length.nc')
    snow_season_length_measure.to_netcdf(savedir_threshold+'measure_snow_season_lenght.nc')


"""apply masks to data"""


snow_doy_racmo = xr.open_dataset(savedir_threshold+'racmo_snow_doy.nc')
snow_doy_measure = xr.open_dataset(savedir_threshold+'measure_snow_doy.nc')

first_masked = xr.DataArray(snow_doy_racmo['first_day'].to_numpy()*land_mask_universal,name='first_day')
last_masked = xr.DataArray(snow_doy_racmo['last_day'].to_numpy()*land_mask_universal,name='last_day')

snow_doy_racmo_masked = xr.merge([first_masked,last_masked])

first_masked = xr.DataArray(snow_doy_measure['first_day'].to_numpy()*land_mask_universal,name='first_day')
last_masked = xr.DataArray(snow_doy_measure['last_day'].to_numpy()*land_mask_universal,name='last_day')

snow_doy_measure_masked = xr.merge([first_masked,last_masked])



"""map with location correspoding to i in areas_int"""

if Location_map == True:


    
    i = 0
    
    
    areas_int = [['Canada',66.11,-70,78],
                 ['Greenland-West',68.32,-52.1,59],
                 ['Greenland-East',71.15,-23.47,59],
                 ['Iceland',65.87,-22.01,60],
                 ['Svalbard',77.99,22.77,83]
                 ]
    
    col , row = calculations.return_index_from_coordinates(areas_int[i][1], areas_int[i][2], racmo_ease)
    
    plt.contourf(land_mask_universal,cmap='Greys')
    plt.colorbar()
    plt.scatter(col,row, color = 'red')
    plt.xlim(250,500)
    plt.ylim(250,500)
    plt.title('Racmo_Land_masked')
    plt.savefig(fig_save_directory+'Map_with_location_'+areas_int[i][0]+'.png')




if Day_of_year_plot == True:
    
    """plot scatter heatmap day of year"""
    
    fig, axs = plt.subplots(1,2,figsize=(10,4),dpi=800)
    
    measure = snow_doy_measure_masked['first_day'][250:500,250:500].values.flatten()
    racmo = snow_doy_racmo_masked['first_day'][250:500,250:500].values.flatten()
    
    racmo = np.nan_to_num(racmo,nan=-1)
    measure = np.nan_to_num(measure,nan=-1)
    
    hist, xedges, yedges = np.histogram2d(measure,racmo,bins=120)
    xidx = np.clip(np.digitize(measure, xedges), 0, hist.shape[0]-1)
    yidx = np.clip(np.digitize(racmo, yedges), 0, hist.shape[1]-1)
    c = hist[xidx, yidx]
    
    axs[0].set_title('First snow-free day in 2001')
    axs[0].scatter(measure, racmo, c=c, s = 1, cmap=plt.cm.RdYlBu_r,norm=mpl.colors.LogNorm())
    axs[0].set_xlabel('Satellite [day of year]')
    axs[0].set_ylabel('Racmo [day of year]')
    axs[0].plot(range(365),range(365),color='black',linestyle=(0,(3,3)),zorder=10,alpha=0.5)
    axs[0].set_xlim(0,365)
    axs[0].set_ylim(0,365)
    axs[0].set_aspect(1)
    axs[0].set_xticks(np.arange(0,365,50))
    axs[0].set_yticks(np.arange(0,365,50))
    
    measure = snow_doy_measure_masked['last_day'][250:500,250:500].values.flatten()
    racmo = snow_doy_racmo_masked['last_day'][250:500,250:500].values.flatten()
    
    racmo = np.nan_to_num(racmo,nan=-1)
    measure = np.nan_to_num(measure,nan=-1)
    
    hist, xedges, yedges = np.histogram2d(measure,racmo,bins=120)
    xidx = np.clip(np.digitize(measure, xedges), 0, hist.shape[0]-1)
    yidx = np.clip(np.digitize(racmo, yedges), 0, hist.shape[1]-1)
    c = hist[xidx, yidx]
    
    
    axs[1].set_title('Last snow-free day in 2001')
    axs[1].scatter(measure, racmo, c=c, s = 1, cmap=plt.cm.RdYlBu_r,norm=mpl.colors.LogNorm())
    axs[1].set_xlabel('Satellite [day of year]')
    axs[1].set_ylabel('Racmo [day of year]')
    axs[1].plot(range(365),range(365),color='black',linestyle=(0,(3,3)),zorder=10,alpha=0.5)
    axs[1].set_xlim(0,365)
    axs[1].set_ylim(0,365)
    axs[1].set_aspect(1)
    axs[1].set_xticks(np.arange(0,365,50))
    axs[1].set_yticks(np.arange(0,365,50))
    
    plt.savefig('fig_save_directory'+'RACMO_MEASURE_DOY_SCATTER.png',dpi=800)
    
    
    """plot only first snow free day"""
    
    vfirst = [0,265]
    vlast= [100,365]
    
    cmap=plt.cm.RdYlBu_r
    levels = 30
    
    fig, axs = plt.subplots(1,2,figsize=(15, 6),sharey=True,sharex=True,dpi=600)
    fig.suptitle('First snow-free day')
    
    axs[0].contourf(snow_doy_measure_masked['first_day'].values,cmap=cmap,vmin=vfirst[0],vmax=vfirst[1],levels=levels)
    axs[0].set_title('Measure')
    axs[0].set_xlim([220,495])
    axs[0].set_ylim([220,495])
    
    first = axs[1].contourf(snow_doy_racmo_masked['first_day'].values,cmap=cmap,vmin=vfirst[0],vmax=vfirst[1],levels=levels)
    axs[1].set_title('Racmo')
    
    plt.tight_layout()
    
    fig.colorbar(first, ax=axs[0:2],location='right',shrink=0.8)
    
    plt.savefig(fig_save_directory+'RACMO_MEASURE_doi_first_snow_map.png',dpi=600)
    
    plt.show()
    
    
    """plot only last snow free day"""
    
    fig, axs = plt.subplots(1,2,figsize=(15, 6),sharey=True,sharex=True,dpi=600)
    fig.suptitle('Last snow-free day')
    
    axs[0].contourf(snow_doy_measure_masked['last_day'].values,cmap=cmap,vmin=vlast[0],vmax=vlast[1],levels=levels)
    axs[0].set_title('Measure')
    axs[0].set_xlim([220,495])
    axs[0].set_ylim([220,495])
    
    last = axs[1].contourf(snow_doy_racmo_masked['last_day'].values,cmap=cmap,vmin=vlast[0],vmax=vlast[1],levels=levels)
    axs[1].set_title('Racmo')
    
    plt.tight_layout()
    
    fig.colorbar(last, ax=axs[0:2],location='right',shrink=0.8)
    
    
    plt.savefig(fig_save_directory+'RACMO_MEASURE_doi_last_snow_map.png',dpi=600)

if Season_length_plot == True: 
    
    """plot season lengths"""
    
    
    
    cols = measure_ease.coords['cols'].values
    rows = measure_ease.coords['rows'].values
        
    snow_season_length_racmo = xr.open_dataset(savedir_threshold+'racmo_snow_season_length.nc')
    snow_season_length_measure = xr.open_dataset(savedir_threshold+'measure_snow_season_lenght.nc')
    
    
    snow_season_length_racmo = snow_season_length_racmo.assign_coords({'latitude':(('rows','cols'),latitudes),'longitude':(('rows','cols'),longitudes),'rows':(('rows'),rows),'cols':(('cols'),cols)})
    snow_season_length_measure = snow_season_length_measure.assign_coords({'latitude':(('rows','cols'),latitudes),'longitude':(('rows','cols'),longitudes),'rows':(('rows'),rows),'cols':(('cols'),cols)})
    
    """Masking"""
    
    snow_season_length_measure['first_season_length_masked'] = snow_season_length_measure['first_season_length'].where(land_mask_universal==1,np.nan)
    snow_season_length_measure['last_season_length_masked'] = snow_season_length_measure['last_season_length'].where(land_mask_universal==1,np.nan)
    
    snow_season_length_racmo['first_season_length_masked'] = snow_season_length_racmo['first_season_length'].where(land_mask_universal==1,np.nan)
    snow_season_length_racmo['last_season_length_masked'] = snow_season_length_racmo['last_season_length'].where(land_mask_universal==1,np.nan)
    
    
    lim = 4000000
    
    vmin = 0
    vmax = 150
    
    cmap='Blues'
    levels = 30
    
    """First Season"""
    
    fig, axs = plt.subplots(1,2, figsize=(12, 5), subplot_kw={'projection':ccrs.NorthPolarStereo()},dpi=800)
    fig.suptitle('First season')
    axs[0].coastlines(resolution='110m',alpha=0.5)
    axs[1].coastlines(resolution='110m',alpha=0.5)
    
    snow_season_length_racmo['first_season_length_masked'].plot(ax=axs[0], transform=ccrs.epsg(6931),cmap=cmap,vmin=vmin,vmax=vmax,levels=levels,add_colorbar=False)
    snow_season_length_measure['first_season_length_masked'].plot(ax=axs[1], transform=ccrs.epsg(6931),cmap=cmap,vmin=vmin,vmax=vmax,levels=levels)
    
    
    ice_mask.plot(ax=axs[0],cmap='Oranges',add_colorbar=False)
    ice_mask.plot(ax=axs[1],cmap='Oranges',add_colorbar=False)
    
    
    axs[0].set_title('Racmo')
    axs[1].set_title('Satellite')
    
    axs[0].set_ylim([-lim,lim])
    axs[0].set_xlim([-lim,lim])
    
    axs[1].set_ylim([-lim,lim])
    axs[1].set_xlim([-lim,lim])
    
    plt.tight_layout()
    
    plt.savefig(fig_save_directory+'RACMO_MEASURE_FIRST_SEASON_LENGTH_MAP.png',dpi=800)
    
    
    """Last Season"""
    
    vmin = 0
    vmax = 100
    
    cmap='Blues'
    levels = 30
    
    fig, axs = plt.subplots(1,2, figsize=(12, 5), subplot_kw={'projection':ccrs.NorthPolarStereo()},dpi=800)
    fig.suptitle('Last season')
    
    axs[0].coastlines(resolution='110m',alpha=0.5)
    axs[1].coastlines(resolution='110m',alpha=0.5)
    
    snow_season_length_racmo['last_season_length_masked'].plot(ax=axs[0], transform=ccrs.epsg(6931),cmap=cmap,vmin=vmin,vmax=vmax,levels=levels,add_colorbar=False)
    snow_season_length_measure['last_season_length_masked'].plot(ax=axs[1], transform=ccrs.epsg(6931),cmap=cmap,vmin=vmin,vmax=vmax,levels=levels)
    
    ice_mask.plot(ax=axs[0],cmap='Oranges',add_colorbar=False)
    ice_mask.plot(ax=axs[1],cmap='Oranges',add_colorbar=False)
    
    axs[0].set_title('Racmo')
    axs[1].set_title('Satellite')
    
    axs[0].set_ylim([-lim,lim])
    axs[0].set_xlim([-lim,lim])
    
    axs[1].set_ylim([-lim,lim])
    axs[1].set_xlim([-lim,lim])
    
    plt.tight_layout()
    
    plt.savefig(fig_save_directory+'RACMO_MEASURE_LAST_SEASON_LENGTH_MAP.png',dpi=800)
    
    
    """Season difference"""
    
    vmin = -50
    vmax = 50
    
    cmap='seismic'
    levels = 30
    
    fig, axs = plt.subplots(1,2, figsize=(12, 5), subplot_kw={'projection':ccrs.NorthPolarStereo()},dpi=800)
    
    axs[0].coastlines(resolution='110m',alpha=0.5)
    axs[1].coastlines(resolution='110m',alpha=0.5)
    
    (snow_season_length_racmo['first_season_length_masked']-snow_season_length_measure['first_season_length_masked']).plot(ax=axs[0], transform=ccrs.epsg(6931),cmap=cmap,vmin=vmin,vmax=vmax,levels=levels,add_colorbar=False)
    (snow_season_length_racmo['last_season_length_masked']-snow_season_length_measure['last_season_length_masked']).plot(ax=axs[1], transform=ccrs.epsg(6931),cmap=cmap,vmin=vmin,vmax=vmax,levels=levels)
    
    ice_mask.plot(ax=axs[0],cmap='Oranges',add_colorbar=False)
    ice_mask.plot(ax=axs[1],cmap='Oranges',add_colorbar=False)
    
    axs[0].set_title('First season')
    axs[1].set_title('Last Season')
    
    axs[0].set_ylim([-lim,lim])
    axs[0].set_xlim([-lim,lim])
    
    axs[1].set_ylim([-lim,lim])
    axs[1].set_xlim([-lim,lim])
    
    plt.tight_layout()
    
    plt.savefig(fig_save_directory+'RACMO_MEASURE_SEASON_LENGTH_DIFFERENCE_MAP.png',dpi=800)
    
    """plot scatter heatmap day of year"""
    

    fig, axs = plt.subplots(1,2,figsize=(10,4),dpi=800)
    
    measure = snow_season_length_measure['first_season_length_masked'][250:500,250:500].values.flatten()
    racmo = snow_season_length_racmo['first_season_length_masked'][250:500,250:500].values.flatten()
    
    racmo = np.nan_to_num(racmo,nan=-1)
    measure = np.nan_to_num(measure,nan=-1)
    
    figlim = 183
    
    hist, xedges, yedges = np.histogram2d(measure,racmo,bins=120)
    xidx = np.clip(np.digitize(measure, xedges), 0, hist.shape[0]-1)
    yidx = np.clip(np.digitize(racmo, yedges), 0, hist.shape[1]-1)
    c = hist[xidx, yidx]
    
    
    RMSE = np.sqrt(np.mean(((racmo-measure)**2)))
    regres =scipy.stats.linregress(racmo,measure)
    
    
    axs[0].set_title('First season length 2001')
    axs[0].scatter(measure, racmo, c=c, s = 1, cmap=plt.cm.RdYlBu_r,norm=mpl.colors.LogNorm())
    axs[0].set_xlabel('Satellite [number of days]')
    axs[0].set_ylabel('Racmo [number of days]')
    axs[0].plot(range(figlim),range(figlim),color='black',linestyle=(0,(3,3)),zorder=10,alpha=0.5)
    axs[0].set_xlim(0,figlim)
    axs[0].set_ylim(0,figlim)
    axs[0].set_aspect(1)
    axs[0].set_xticks(np.arange(0,figlim,50))
    axs[0].set_yticks(np.arange(0,figlim,50))
    axs[0].annotate(('RMSE:'+str(np.round(RMSE,1))),xy=(100,50))
    axs[0].annotate(('Slope:'+str(np.round(regres.slope,3))),xy=(100,40))  
    axs[0].annotate(('CC:'+str(np.round(regres.rvalue,3))),xy=(100,30))
    
    measure = snow_season_length_measure['last_season_length_masked'][250:500,250:500].values.flatten()
    racmo = snow_season_length_racmo['last_season_length_masked'][250:500,250:500].values.flatten()
    
    racmo = np.nan_to_num(racmo,nan=-1)
    measure = np.nan_to_num(measure,nan=-1)
    
    RMSE = np.sqrt(np.mean(((racmo-measure)**2)))
    regres =scipy.stats.linregress(racmo,measure)
    
    hist, xedges, yedges = np.histogram2d(measure,racmo,bins=120)
    xidx = np.clip(np.digitize(measure, xedges), 0, hist.shape[0]-1)
    yidx = np.clip(np.digitize(racmo, yedges), 0, hist.shape[1]-1)
    c = hist[xidx, yidx]
    
    
    axs[1].set_title('Last season length 2001')
    axs[1].scatter(measure, racmo, c=c, s = 1, cmap=plt.cm.RdYlBu_r,norm=mpl.colors.LogNorm())
    axs[1].set_xlabel('Satellite [number of days]')
    axs[1].set_ylabel('Racmo [number of days]')
    axs[1].plot(range(figlim),range(figlim),color='black',linestyle=(0,(3,3)),zorder=10,alpha=0.5)
    axs[1].set_xlim(0,figlim)
    axs[1].set_ylim(0,figlim)
    axs[1].set_aspect(1)
    axs[1].set_xticks(np.arange(0,figlim,50))
    axs[1].set_yticks(np.arange(0,figlim,50))
    axs[1].annotate(('RMSE:'+str(np.round(RMSE,1))),xy=(100,50))
    axs[1].annotate(('Slope:'+str(np.round(regres.slope,3))),xy=(100,40))  
    axs[1].annotate(('CC:'+str(np.round(regres.rvalue,3))),xy=(100,30))
    
    plt.savefig(fig_save_directory+'RACMO_MEASURE_SEASON_LENGTH_SCATTER.png',dpi=800)