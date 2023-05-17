#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 11:39:56 2023

@author: tijmen
"""

import os

import pandas as pd
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

years = [2005, 2006, 2007]

for year in years:
    year = str(year)

    print('making plots for: ' + year)

    t_start = year+'-01-01'
    t_stop = year+'-12-31'


    racmo_arctic_data_directory = datadir+'RACMO_2.4/PXARC11/2001_new/'
    directory4_DEFAULT = racmo_arctic_data_directory + 'NC_DEFAULT/'
    directory4_MD = racmo_arctic_data_directory + 'NC_MD/'
    remapdir = datadir+'remap/'
    modis_data_directory = datadir+'MODIS/'
    aws_directory = datadir + 'Aws_data/'
    in_situ_data_directory_year_calculated = datadir + 'In_situ_data/' + year + '/Calculated/'

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

    MODIS = xr.open_dataset(modis_data_directory+year+'_RCG.nc')

    def monthly_stationdata_import(directory):
        yearly_df = pd.DataFrame()

        for i in range(12):
            path = (directory+'month_'+str(i+1)+'/stationdata.csv')
            monthly_df = pd.read_csv(path, index_col=0)
            yearly_df = pd.concat([yearly_df, monthly_df])

        return(yearly_df)

    IN_SITU_TAS = monthly_stationdata_import('/Volumes/Tijmen/Master-Thesis/Data/In_situ_data/'+year+'/Calculated/air_temperature/')
    IN_SITU_SNWD = monthly_stationdata_import('/Volumes/Tijmen/Master-Thesis/Data/In_situ_data/'+year+'/Calculated/snow_depth/')
    IN_SITU_PREC = monthly_stationdata_import('/Volumes/Tijmen/Master-Thesis/Data/In_situ_data/'+year+'/Calculated/accumulated_precipitation/')

    #%%

    save_directory = '/Users/tijmen/Desktop/Figures_Thesis/'

    areas_int = [['SODANKYLA AWS GSN',67.368, 26.633, 0],
                 ['NORRBACK', 64.71, 17.72, 0]
                 ]



    i = 0
    AWS = pd.read_csv(aws_directory+'SODANKYLA.csv')
    AWS.rename(columns={'m': 'month', 'd': 'day'}, inplace=True)
    AWS['datetime'] = pd.to_datetime(AWS[['Year', 'month', 'day']])

    # Set the datetime column as the index of the DataFrame
    AWS = AWS.set_index('datetime')


    index = calculations.return_index_from_coordinates(areas_int[i][1],areas_int[i][2],data.Variable_Import(directory4_MD, 'rlds'))
    modis_index =  calculations.return_index_from_coordinates(areas_int[i][1],areas_int[i][2], MODIS)

    modis_rlat = modis_index[1]
    modis_rlon = modis_index[0]
    rlat = index[1]
    rlon = index[0]

    fig, axs = plt.subplots(2,3, figsize=(40, 20))
    fig.suptitle(areas_int[i][0],fontsize=20)

    """RACMO 4 PLOTTING"""

    axs[0,0].plot(RC4_LW_D.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='orange',label='Down')
    axs[0,0].plot(RC4_LW_U.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='blue',label='Up')
    axs[0,0].plot((RC4_LW_U+RC4_LW_D).isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='black',label='net')
    axs[0,0].set_title('Longwave')



    axs[1,0].plot(RC4_SW_D.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start, t_stop)).values,color='orange',label='Down')
    axs[1,0].plot(RC4_SW_U.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start, t_stop)).values,color='blue',label='Up')
    axs[1,0].plot((RC4_SW_U+RC4_SW_D).isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='black',label='net')
    axs[1,0].set_title('Shortwave')


    axs[0,1].plot(RC4_SNWD.isel(rlon=rlon, rlat=rlat).sel(time=slice(t_start, t_stop)).values*100, color='Blue', label='Racmo')
    axs[0,1].plot(IN_SITU_SNWD[areas_int[i][0]], color='red', label='In_situ')
    axs[0,1].set_title('Snowheight')

    axs[0,2].plot(IN_SITU_PREC[areas_int[i][0]], color='red', label='In_situ')
    axs[0,2].set_title('Precipitation')

    axs[1,1].plot((abs(RC4_SW_U)/abs(RC4_SW_D)).isel(rlon=rlon, rlat=rlat).sel(time=slice(t_start, t_stop)).values, 'Blue', linestyle=(0, (4, 4)), label='Racmo SW_U/SW_D')
    axs[1,1].plot(MODIS['Albedo'].isel(rlon=modis_rlon, rlat=modis_rlat).sel(time=slice(t_start, t_stop)).values, color='Green', label='Modis')
    axs[1,1].set_title('Albedo')

    axs[1,2].plot(RC4_TAS.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).values,color='Blue', label='Racmo')
    axs[1,2].plot(IN_SITU_TAS[areas_int[i][0]], color='red', label='In_situ')
    axs[1,2].set_title('Temperature')


    """apply to all pannels"""

    for y in range(len(axs[0,:])):
        for x in range(len(axs[:,0])):
            ylims = axs[x,y].get_ylim()
            axs[x,y].vlines(areas_int[i][3], ylims[0],ylims[1],color='black',linestyle=(0,(3,3)))
            axs[x,y].legend()

    os.makedirs(save_directory+'In_situ_TS/'+year+'/', exist_ok=True)

    plt.savefig(save_directory+'In_situ_TS/'+year+'/'+areas_int[i][0]+'_TS_plot')


    if False == True:
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.coastlines(alpha=0.5,resolution='50m')


        for i in range(len(areas_int[0])):
            ax.scatter(areas_int[i][2],areas_int[i][1],zorder=10,s=30,label=areas_int[i][0])


        plt.xlim((-180,180))
        plt.ylim((0,90))
        plt.legend()
        plt.show()

        #%%
