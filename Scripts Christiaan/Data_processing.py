#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 15:02:34 2023

@author: tijmen
"""

desktop = False
laptop = True




import os
import numpy as np
import pandas as pd
import datetime
import xarray as xr
import statistics
from shapely.geometry import Polygon, Point
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import geopandas as gpd
import math

if laptop:
    os.chdir('/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Python_scripts')

if desktop:
    os.chdir('H:\Documenten\Master\Master_thesis\Python_scripts')

import Thesis_Functions.calculations as Calculations
import Thesis_Functions.data as Data


"""Variables"""

years = [2002, 2003, 2004, 2005]            #list of years where data should be proccessed over, entire year is processed. Data should exist in format as specified
months = [1,2,3,4,5,6,7,8,9,10,11,12]                       #list of months to process
breakdate = '-07-01'                    #Split date between accumulation and melt season
days_missing_limit = 5                  #Maximum number of missing days before station is discarted (MAKE MORE REFINED FILTER)
tilefrac = 'tilefrac5'                  #tilefraction over which calculations are done
tilefrac_threshold = 0                  #threshold for tilefraction
moving_average_window = 3
constant_date = '-04-30'
"""Calculation control"""

select_stations = False                 # Select stations that are in arctic domain
tilefrac_select = False                 # select stations whithin a certain tilefraction
calc_stationdata = False                # Extract station snowdepth data and save to csv per station
interpolate_stationdata = False         # Fill missing values
fill_nan = False                        # Interpolate
monthly_data = False                    # Extract montly data and save to directories according to structure: /Year/Month/variable.nc
monthly_data_racmo_only = False         # extract monthly data at in_situ measurement locations in racmo only
select_stations_area = False            # select stations that are in focus area
monthly_statistics = False              # get monthly statistics (bias, mean, RMSE)
racmo_snowextent = False                # get total snow covered tilefraction
combine_snow_extend = False             # Combine snow extend netcdfs in one file
snow_extend_statistics = False          # calculate length of melt and accumulation season
albedo_extraction = False               # extract modis albedo at in_situ measurement locations
mask_albedo = False                     # select only tiles that are completely covered with snow
mask_albedo_constant = True
aws_data_modification = False           # AWS data from finnland saved as daily mean
aws_weather = False                     # AWS weather data to yearly files
racmo_clear_sky = False                 # Calculate racmo clear sky albedo from clear sky radiative fluxes
identify_albedo_events_aws = False      # Identify albedo events in racmo timeseries

"""File names"""

measure_filename='/Measure_merged.nc' #Filename for combined measure dataset

"""Variable control"""

Snowdepth = True
Surface_temp = False
Precipitation = False
Albedo = False

"""Directories"""

if laptop:
    print('Directory structure: laptop')
    datadir = "/Volumes/Tijmen/Master-Thesis/Data/"
if desktop:
    print('Directory structure: desktop')
    datadir = "E:/Master-Thesis/Data/"

racmo_v1 = False
racmo_v2 = True

if racmo_v1:

    in_situ_data_directory = datadir+'In_situ_data/'
    racmo_arctic_data_directory = datadir+'RACMO_2.4/PXARC11/2001_new/'
    snow_cover_extend_measure_dir = datadir+'Snow_cover_Measure/'
    mask_directory = datadir+'Mask/'
    remapdir = datadir+'Remap/'
    snow_cover_analysis_dir = datadir+'Snow_cover_analyses/Snow_cover_ease_grid/'
    download_measure_dir = 'Download_3-4/'
    modis_data_directory = datadir+'MODIS/'
    aws_directory = datadir + 'Aws_data/'

    racmo_filename_additive = '.KNMI-2001.PXARC11.RACMO24_1_complete6_UAR_q_noice_khalo6_era5q.DD.nc'

if racmo_v2:

    in_situ_data_directory = datadir+'In_situ_data/'
    racmo_arctic_data_directory = datadir+'RACMO_2.4/PXARC11/2001_v2/'
    snow_cover_extend_measure_dir = datadir+'Snow_cover_Measure/'
    mask_directory = datadir+'Mask/'
    remapdir = datadir+'Remap/'
    snow_cover_analysis_dir = datadir+'Snow_cover_analyses/Snow_cover_ease_grid/'
    download_measure_dir = 'Download_3-4/'
    modis_data_directory = datadir+'MODIS/'
    aws_directory = datadir + 'Aws_data/'

    racmo_filename_additive = '.KNMI-2001.PXARC11.RACMO24_1_complete6_UAR_q_era5q_tijmen.DD.nc'


"""Bounding boxes for location extraction"""

polygon_norway = Polygon([(4.0, 55.0), (4.0, 65.0), (31.0, 75.0), (31.0, 70.0), (17.0, 66.0), (12.0, 58.0)])

polygon_flat_europe = Polygon([(4.0, 55.0), (12.0, 58.0), (17.0, 66.0), (31.0, 70.0), (31.0, 55.0)])

polygon_syberia = Polygon([(31.0, 55.0), (31.0, 75.0), (180.0, 75.0), (190.0, 66.0), (180.0, 55.0)])

polygon_alaska = Polygon([(-170.0, 55.0), (-170.0, 65.0), (-160.0, 75.0), (-160.0, 75.0), (-140.0, 55.0)])

polygon_canada = Polygon([(-130.0, 65.0), (-130.0, 85.0), (-60.0, 85.0), (-80.0, 75.0), (-60.0, 65.0)])

"""Start of yearly calculations:"""

if Snowdepth:
    print('Calculation variable: Snowdepth')

    in_situ_variable = 'snow_depth'
    racmo_variable = 'sndp'
    racmo_filename = 'NC_DEFAULT/'+racmo_variable+racmo_filename_additive


if Surface_temp:
    print('Calculation variable: Surface temperature')

    in_situ_variable = 'air_temperature'
    racmo_variable = 'tas'
    racmo_filename = 'NC_DEFAULT/'+racmo_variable+racmo_filename_additive


if Precipitation:

    in_situ_variable = 'accumulated_precipitation'
    racmo_variable = 'pr'
    racmo_filename = 'NC_DEFAULT/'+racmo_variable+racmo_filename_additive

if Albedo:
    in_situ_variable = '-'
    racmo_variable = 'Clearsky_albedo_calculated'
    racmo_filename = 'NC_MC/'+racmo_variable+'.nc'


if tilefrac_select:

    tilefrac_filename = 'NC_DEFAULT/'+tilefrac+racmo_filename_additive

for _ , year in enumerate(years):

    year = str(year)
    print('Starting calculations for year: '+ year)
    """Year specific directories"""
    
    in_situ_data_directory_year = in_situ_data_directory+year+'/'
    in_situ_data_directory_year_calculated = in_situ_data_directory+year+'/Calculated/'

    """If directory does not yet exists, make directories"""
    
    os.makedirs(in_situ_data_directory_year_calculated, exist_ok=True)

    """Select stations that are in arctic domain"""

    if select_stations:

        racmo_24_arc_snowheight = xr.open_dataset(racmo_arctic_data_directory + racmo_filename).sel(time=slice(year + '-01-01', year + '-12-31'))[racmo_variable]
        rlats = racmo_24_arc_snowheight['rlat'].values
        rlons = racmo_24_arc_snowheight['rlon'].values

        in_situ_data_directory_yearlist = os.listdir(in_situ_data_directory_year+'/Raw_data/')
        openfiles = []

        for i, v in enumerate(in_situ_data_directory_yearlist):
            try:
                in_situ_data_directory_yearlist[i].split('.')[1] == 'csv'
            except:
                None
            else:
                if in_situ_data_directory_yearlist[i].split('.')[1] == 'csv':
                    openfiles.append(v)

        observations = []

        for i, v in enumerate(openfiles):
            observations.append(pd.read_csv(in_situ_data_directory_year+'/Raw_data/' + v))

        observation = pd.concat(observations)

        stations = np.sort(list(set(observation['station_name'])))
        indexes = ['latitude', 'longitude', 'rlat', 'rlon', 'in_arctic_domain']
        station_stats = pd.DataFrame(index=indexes, columns=stations)

        print('Selecting stations that are in racmo range')

        for i, v in enumerate(stations):

            station_data = observation.loc[observation['station_name']==v].iloc[0]

            lat = station_data['latitude']
            lon = station_data['longitude']

            x, y = Calculations.return_index_from_coordinates(lat, lon, racmo_24_arc_snowheight)

            rlat = rlats[y]
            rlon = rlons[x]

            if x != 0 and x != len(racmo_24_arc_snowheight['rlon']) - 1 and y != 0 and y != len(racmo_24_arc_snowheight['rlat']) - 1:
                station_stats.loc['in_arctic_domain',v]=True
                station_stats.loc['latitude', v] = lat
                station_stats.loc['longitude', v] = lon
                station_stats.loc['rlat', v] = rlat
                station_stats.loc['rlon', v] = rlon
            else:
                station_stats.drop(columns=v,inplace=True)

        station_stats.to_csv(in_situ_data_directory_year_calculated + 'station_in_arctic_domain_' + year + '.csv')

    """Select stations within a certain tile with feature 'tilefrac' in it"""
    if tilefrac_select:

        racmo_tilefrac = xr.open_dataset(racmo_arctic_data_directory + tilefrac_filename).sel(time=slice(year + '-01-01', year + '-12-31'))[tilefrac]

        station_stats = pd.read_csv(in_situ_data_directory_year_calculated + 'station_in_arctic_domain_' + year + '.csv', index_col=0)

        print('Selecting stations that are in '+tilefrac)

        for i, v in enumerate(station_stats.columns):

            rlat = station_stats.loc['rlat'][v]
            rlon = station_stats.loc['rlon'][v]

            if racmo_tilefrac.sel(rlon=rlon, rlat=rlat).max('time').values[0] > tilefrac_threshold:

                station_stats.drop(columns=v,inplace=True)

        print(str(len(station_stats.columns)) + ' stations found in ' + tilefrac)

        station_stats.to_csv(in_situ_data_directory_year_calculated + 'station_'+tilefrac + '_' + year + '.csv')

    """Reformat raw data into station specific data, only when true"""
    
    if calc_stationdata:
        
        print('Transforming stationdata: per station yearly record')
        
        in_situ_data_directory_yearlist = os.listdir(in_situ_data_directory_year+'/Raw_data/')
        openfiles = []
        
        for i,v in enumerate(in_situ_data_directory_yearlist):
            try: in_situ_data_directory_yearlist[i].split('.')[1] == 'csv'
            
            except:
                None
            else: 
                if in_situ_data_directory_yearlist[i].split('.')[1] == 'csv':
    
                    openfiles.append(v)
        
        observations = []
        
        for i , v in enumerate(openfiles):
            
            observations.append(pd.read_csv(in_situ_data_directory_year+'/Raw_data/'+v))

        observation = pd.concat(observations)

        station_stats = pd.read_csv(in_situ_data_directory_year_calculated + 'station_in_arctic_domain_' + year + '.csv',index_col=0)

        stations = station_stats.columns

        dateindex = pd.date_range(year+"-01-01", year+"-12-31", freq="D")

        stationarray = pd.DataFrame(index=dateindex, columns =stations)

        for i,v in enumerate(stations):
            
            station_snow_depth = observation.loc[(observation['station_name'] == v) & (observation['observed_variable'] == in_situ_variable)]
            for i2 , v2 in enumerate(station_snow_depth['observation_value']):
                try: stationarray.loc[station_snow_depth['date_time'].iloc[i2].split(' ')[0],v] = v2
                
                except: stationarray.loc[station_snow_depth['date_time'].iloc[i2],v] = v2
        
        stationarray.to_csv(in_situ_data_directory_year_calculated+'stations_daily_'+in_situ_variable+'_'+ year + '.csv')

    """Interpolate and select data with a bigger gap then days_missing_limit"""
    
    if interpolate_stationdata:
        
        print('Interpolating stationdata')
       
        stationarray = pd.read_csv(in_situ_data_directory_year_calculated+'stations_daily_'+in_situ_variable+'_'+ year + '.csv',index_col=0)
        stations = np.array(stationarray.columns)
        
        for i,v in enumerate(stations):
    
            stationarray[v].interpolate(method='linear', inplace=True, limit = days_missing_limit)
    
        stationarray.to_csv(in_situ_data_directory_year_calculated+'stations_daily_'+in_situ_variable+'_interpolated_'+year+'.csv')

    """Fill nan's in in situ snowheight data in summer"""

    if fill_nan:

        print("Filling nans:")

        """In situ data"""

        stationarray = pd.read_csv(
            in_situ_data_directory_year_calculated + 'stations_daily_'+in_situ_variable+'_' + year + '.csv', index_col=0)

        station_stats = pd.read_csv(
            in_situ_data_directory_year_calculated + 'station_in_arctic_domain_' + year + '.csv', index_col=0)

        stationarray.index = pd.to_datetime(stationarray.index)

        # Create a new column to identify summer months
        stationarray['is_summer'] = ((stationarray.index.month >= 5) & (stationarray.index.month <= 9))

        for station in stationarray.columns:
            stationarray.loc[(stationarray.index.month >= 5) & (stationarray.index.month <= 9), station] = \
            stationarray[station][stationarray['is_summer']].apply(lambda x: 0 if math.isnan(x) else x)
        stationarray = stationarray.drop('is_summer', axis=1)

        stationarray.to_csv(
            in_situ_data_directory_year_calculated + 'stations_daily_'+in_situ_variable+'_no_nan_' + year + '.csv')

    """Station specific data to csv"""

    if monthly_data:
        
        print("Test: saving racmo and in-situ data as monthly files, month:")
        
        """In situ data"""
        
        stationarray = pd.read_csv(in_situ_data_directory_year_calculated+'stations_daily_'+in_situ_variable+'_'+year+'.csv', index_col=0)
        racmo_24_arc_snowheight = xr.open_dataset(racmo_arctic_data_directory+racmo_filename) #.sel(time=slice(year +'-01-01',year+'-12-31'))[racmo_variable]
        station_stats = pd.read_csv(in_situ_data_directory_year_calculated+'station_in_arctic_domain_'+year+'.csv', index_col=0)
        
        rlats = station_stats.loc['rlat']
        rlons = station_stats.loc['rlon']

        for month in months:

            if len(racmo_24_arc_snowheight.time.values) < 364:
                print('year '+year+' is not completely included in the racmo dataset, please select different dataset.')
                break
            
            print('Month '+ str(month))

            start_date = pd.to_datetime(f'{year}-{month}-01')
            end_date = pd.to_datetime(f'{year}-{month}-{pd.Period(start_date, freq="M").days_in_month}')

            """get data for given month"""

            month_racmo = racmo_24_arc_snowheight.sel(time=slice(start_date, end_date))
            month_in_situ = stationarray.loc[str(start_date):str(end_date)]

            if len(month_racmo.time) != len(month_in_situ.index):
                month_in_situ = stationarray.loc[str(start_date-datetime.timedelta(days=1)):str(end_date)]

            if len(month_racmo.time) != len(month_in_situ.index):
                print('Length of racmo and in situ dates is not the same')
                break

            month_racmo_per_station = pd.DataFrame(index=month_in_situ.index, columns=station_stats.columns)

            for i, v in enumerate(station_stats.columns):

                month_racmo_per_station[v] = month_racmo.sel(rlat=station_stats.loc['rlat',v],rlon=station_stats.loc['rlon',v])


            """Check for directory to exist"""
            
            monthdir_in_situ = in_situ_data_directory_year_calculated +'/'+in_situ_variable+'/month_'+str(month)
            os.makedirs(monthdir_in_situ,exist_ok=True)
        
            monthdir_racmo = racmo_arctic_data_directory+racmo_variable+'/'+year+'/month_'+str(month)
            os.makedirs(monthdir_racmo,exist_ok=True)
            
            """Get racmo snowheight for same locations"""

            month_in_situ.to_csv(monthdir_in_situ+'/stationdata.csv')
            month_racmo_per_station.to_csv(monthdir_racmo+'/stationdata.csv')

    """Select stations in a certain latitude-longitude box"""
    if select_stations_area:

        print('Selecting stations in polygon')
        # Load the table containing the latitude and longitude information
        station_stats = pd.read_csv(in_situ_data_directory_year_calculated+'station_statistics_'+year+'.csv', index_col=0)

        station_arctic_domain = pd.read_csv(
            in_situ_data_directory_year_calculated + 'station_in_arctic_domain_' + year + '.csv',
            index_col=0)

        station_stats = station_stats[station_arctic_domain.columns.values]

        """Norway"""

        points = gpd.GeoDataFrame(
            geometry=gpd.points_from_xy(station_stats.loc['longitude'], station_stats.loc['latitude']))
        norway_points = points.within(polygon_norway)
        stations = station_stats.columns
        station_stats[stations[norway_points]].to_csv(
            in_situ_data_directory_year_calculated + 'stations_in_norway_' + year + '.csv')

        """Syberia"""

        syberia_points = points.within(polygon_syberia)
        stations = station_stats.columns
        station_stats[stations[syberia_points]].to_csv(
            in_situ_data_directory_year_calculated + 'stations_in_syberia_' + year + '.csv')

        """flat_europe"""

        flat_europe_points = points.within(polygon_flat_europe)
        stations = station_stats.columns
        station_stats[stations[flat_europe_points]].to_csv(
            in_situ_data_directory_year_calculated + 'stations_in_flat_europe_' + year + '.csv')

        """alaska"""

        alaska_points = points.within(polygon_alaska)
        stations = station_stats.columns
        station_stats[stations[alaska_points]].to_csv(
            in_situ_data_directory_year_calculated + 'stations_in_alaska_' + year + '.csv')

        """canada"""


        canada_points = points.within(polygon_canada)
        stations = station_stats.columns
        station_stats[stations[canada_points]].to_csv(
            in_situ_data_directory_year_calculated + 'stations_in_canada_' + year + '.csv')

    """Calculate statistics for each monthly file"""

    if monthly_statistics:

        print('calculating monthly statistics:')

        station_stats = pd.read_csv(in_situ_data_directory_year_calculated+'station_in_arctic_domain_'+year+'.csv', index_col=0)

        stations = station_stats.columns
        indexes = ['racmo_mean','in_situ_mean','variance_racmo','variance_in_situ','bias']



        for month in months:

            print('Month ' + str(month))
            statistics_stations = pd.DataFrame(index=indexes, columns=stations)

            monthdir_in_situ = in_situ_data_directory_year_calculated + in_situ_variable + '/month_' + str(month)
            monthdir_racmo = racmo_arctic_data_directory+racmo_variable+'/'+year+'/month_'+str(month)

            in_situ = pd.read_csv(monthdir_in_situ + '/stationdata.csv',index_col=0)
            racmo = pd.read_csv(monthdir_racmo + '/stationdata.csv',index_col=0)

            racmo = racmo.mul(100)

            statistics_stations.loc['racmo_mean'] = racmo[stations].mean(skipna=True)
            statistics_stations.loc['in_situ_mean'] = in_situ[stations].mean(skipna=True)

            statistics_stations.loc['variance_racmo'] = racmo[stations].var()
            statistics_stations.loc['variance_in_situ'] = in_situ[stations].var()

            statistics_stations.loc['bias'] = racmo[stations].subtract(in_situ[stations]).mean()

            statistics_stations.to_csv(monthdir_in_situ+'/Calculated_statistics.csv')

            del statistics_stations


    """Get snowextend from racmo data [Boolian]"""

    if racmo_snowextent:

        t_start = year + '-01-01'
        t_stop = year + '-12-30'

        tilefrac5 = xr.open_dataset(racmo_arctic_data_directory+'NC_DEFAULT/tilefrac5'+racmo_filename_additive)['tilefrac5'].sel(time=slice(t_start, t_stop))
        tilefrac7 = xr.open_dataset(racmo_arctic_data_directory+'NC_DEFAULT/tilefrac7'+racmo_filename_additive)['tilefrac7'].sel(time=slice(t_start, t_stop))

        tilefrac_57 = (tilefrac5 + tilefrac7)
        tilefrac_57 = tilefrac_57.rename('Snowextend Racmo')

        tilefrac_57.to_netcdf(remapdir + 'RACMO2.4_V2_Snowextend_' + year + '_RP_grid.nc')

    """Combine snowextend daily files to yearly files"""

    if combine_snow_extend:

        directories = os.listdir(snow_cover_extend_measure_dir + download_measure_dir)
        directories = sorted((f for f in directories if not f.startswith(".")), key=str.lower)

        data = []

        for directory in directories:
            files = os.listdir(snow_cover_extend_measure_dir + download_measure_dir + directory)
            file = sorted(file for file in files if not file.startswith(".") and file.endswith(".nc"))

            data.append(
                xr.open_dataset(snow_cover_extend_measure_dir + download_measure_dir + directory + '/' + file[0]))

        'merge days to one file'

        datamerge = xr.concat(data, 'time')

        'change bytes data to snow [1] or no snow [0]'

        datamerge['merged_snow_cover_extent'] = datamerge['merged_snow_cover_extent'].where(
            datamerge['merged_snow_cover_extent'] < 14, np.nan) * 0 + 1
        datamerge['merged_snow_cover_extent'] = datamerge['merged_snow_cover_extent'].where(
            datamerge['merged_snow_cover_extent'] == 1, 0)

        datamerge['modis_cloud_gap_filled_snow_cover_extent'] = datamerge[
                                                                    'modis_cloud_gap_filled_snow_cover_extent'].where(
            datamerge['modis_cloud_gap_filled_snow_cover_extent'] == 10, np.nan) * 0 + 1
        datamerge['modis_cloud_gap_filled_snow_cover_extent'] = datamerge[
            'modis_cloud_gap_filled_snow_cover_extent'].where(
            datamerge['modis_cloud_gap_filled_snow_cover_extent'] == 1, 0)

        t_start = year + '-01-01'
        t_stop = year + '-12-30'

        datamerge = datamerge.sel(time=slice(t_start, t_stop))

        datamerge.to_netcdf(snow_cover_extend_measure_dir + year + measure_filename)

    if snow_extend_statistics:

        print('Calculating snow_extent statistics')

        snow_cover_measure = xr.open_dataset(snow_cover_analysis_dir + year + '/Measure.nc')['merged_snow_cover_extent']
        snow_cover_racmo = xr.open_dataset(snow_cover_analysis_dir + year + '/RACMO_V2.nc')['Snowextend Racmo']

        snow_cover_racmo = snow_cover_racmo.where(snow_cover_racmo > 0.4, 0)
        snow_cover_racmo = snow_cover_racmo.where(snow_cover_racmo == 0, 1)

        snow_cover_measure_melt = snow_cover_measure.sel(time=slice(year+'-01-01', year+breakdate)).sum(dim='time')
        snow_cover_racmo_melt = snow_cover_racmo.sel(time=slice(year+'-01-01', year+breakdate)).sum(dim='time')

        snow_cover_measure_acc = snow_cover_measure.sel(time=slice(year+breakdate, year+'-12-31',)).sum(dim='time')
        snow_cover_racmo_acc = snow_cover_racmo.sel(time=slice(year+breakdate, year+'-12-31',)).sum(dim='time')

        snow_cover_measure_melt.to_netcdf(snow_cover_analysis_dir + year + '/measure_melt_season.nc')
        snow_cover_racmo_melt.to_netcdf(snow_cover_analysis_dir + year + '/racmo_v2_melt_season.nc')

        snow_cover_measure_acc.to_netcdf(snow_cover_analysis_dir + year + '/measure_acc_season.nc')
        snow_cover_racmo_acc.to_netcdf(snow_cover_analysis_dir + year + '/racmo_v2_acc_season.nc')

    if monthly_data_racmo_only:

        print("Test: saving racmo monthly files, month:")

        racmo_24_arc_variable = xr.open_dataset(racmo_arctic_data_directory + racmo_filename)
        station_stats = pd.read_csv(
            in_situ_data_directory_year_calculated + 'station_in_arctic_domain_' + year + '.csv', index_col=0)

        rlats = station_stats.loc['rlat']
        rlons = station_stats.loc['rlon']

        for month in months:

            if len(racmo_24_arc_variable.time.values) < 364:
                print(
                    'year ' + year + ' is not completely included in the racmo dataset, please select different dataset.')
                break

            print('Month ' + str(month))

            start_date = pd.to_datetime(f'{year}-{month}-01')
            end_date = pd.to_datetime(f'{year}-{month}-{pd.Period(start_date, freq="M").days_in_month}')

            """get data for given month"""
            month_in_situ = pd.read_csv(in_situ_data_directory_year_calculated+'snow_depth/month_'+str(month)+'/stationdata.csv', index_col=0)

            month_racmo = racmo_24_arc_variable.sel(time=slice(start_date, end_date)).squeeze()[racmo_variable]
            month_racmo_per_station = pd.DataFrame(index=month_in_situ.index, columns=station_stats.columns)

            for i, v in enumerate(station_stats.columns):
                month_racmo_per_station[v] = month_racmo.sel(rlat=station_stats.loc['rlat', v],
                                                             rlon=station_stats.loc['rlon', v])

            """Check for directory to exist"""

            monthdir_racmo = racmo_arctic_data_directory + racmo_variable + '/' + year + '/month_' + str(month)
            os.makedirs(monthdir_racmo, exist_ok=True)
            month_racmo_per_station.to_csv(monthdir_racmo + '/stationdata.csv')

    if albedo_extraction:

        print('Extracting modis albedo data at station locations per month')

        modis = xr.open_dataset(modis_data_directory+year+'_RCG.nc')
        station_stats = pd.read_csv(
            in_situ_data_directory_year_calculated + 'station_in_arctic_domain_' + year + '.csv', index_col=0)

        for month in months:

            if len(modis.time.values) < 364:
                print(
                    'year ' + year + ' is not completely included in the modis dataset, please select different dataset.')
                break

            print('Month ' + str(month))

            start_date = pd.to_datetime(f'{year}-{month}-01')
            end_date = pd.to_datetime(f'{year}-{month}-{pd.Period(start_date, freq="M").days_in_month}')

            """get data for given month"""
            month_in_situ = pd.read_csv(
                in_situ_data_directory_year_calculated + 'snow_depth/month_' + str(month) + '/stationdata.csv',
                index_col=0)

            month_modis = modis.sel(time=slice(start_date, end_date))
            month_modis_per_station = pd.DataFrame(index=month_in_situ.index, columns=station_stats.columns)

            for i, v in enumerate(station_stats.columns):
                coord = [station_stats.loc['latitude', v], station_stats.loc['longitude', v]]
                coordarray = (abs(month_modis['lon'].astype(float)-float(coord[1])) + abs(month_modis['lat'].astype(float)-float(coord[0])))
                modis_index = np.where(coordarray == np.min(coordarray))

                month_modis_per_station[v] = month_modis['Albedo'].isel(rlat=modis_index[0],
                                                                        rlon=modis_index[1]).squeeze().values

            """Check for directory to exist"""

            monthdir_modis = modis_data_directory + year + '/month_' + str(month)
            os.makedirs(monthdir_modis, exist_ok=True)
            month_modis_per_station.to_csv(monthdir_modis + '/stationdata.csv')

    if mask_albedo:

        print('Masking albedo for full snow cover tiles')

        tilefrac5 = xr.open_dataset(racmo_arctic_data_directory + 'NC_DEFAULT/tilefrac5'+racmo_filename_additive).sel(
            time=slice(year + '-01-01', year + '-12-31'))['tilefrac5']
        tilefrac7 = xr.open_dataset(racmo_arctic_data_directory + 'NC_DEFAULT/tilefrac7'+racmo_filename_additive).sel(
            time=slice(year + '-01-01', year + '-12-31'))['tilefrac7']
        zsa = xr.open_dataset(racmo_arctic_data_directory + 'NC_DEFAULT/cosznt'+racmo_filename_additive).sel(
            time=slice(year + '-01-01', year + '-12-31'))['cosznt']

        snowcover = tilefrac5 + tilefrac7

        del(tilefrac5)
        del(tilefrac7)

        modis_albedo = xr.open_dataset(modis_data_directory + year+'_RCG.nc')['Albedo']

        modis_albedo_no = modis_albedo.where(modis_albedo['lat'] < 75)
        modis_albedo_no = modis_albedo_no.where(modis_albedo_no.values > 0.05)
        modis_albedo_no = modis_albedo_no.where(zsa.squeeze().values < np.cos(np.radians(55)))
        modis_albedo_rolmean = modis_albedo_no.rolling(time=moving_average_window, center=True).mean()
        modis_albedo_masked = modis_albedo_rolmean.where(snowcover.squeeze().values == 1)
        del(modis_albedo_rolmean)
        modis_albedo_masked.to_netcdf(modis_data_directory + year + '_RCG_masked_new.nc')
        del(modis_albedo_masked)

        racmo_albedo = xr.open_dataset(racmo_arctic_data_directory+'NC_MD/Clearsky_albedo_calculated.nc').squeeze().sel(
            time=slice(year + '-01-01', year + '-12-31'))['Clear-sky_albedo']
        racmo_albedo = racmo_albedo.where(modis_albedo.values > 0.05, np.nan)
        racmo_albedo = racmo_albedo.where(racmo_albedo['lat'] <75)
        racmo_albedo = racmo_albedo.where(zsa.squeeze().values < np.cos(np.radians(55)))
        del(modis_albedo)
        racmo_albedo_rolmean = racmo_albedo.rolling(time=moving_average_window, center=True).mean()
        del(racmo_albedo)
        racmo_albedo_masked = racmo_albedo_rolmean.where(snowcover == 1)
        del (racmo_albedo_rolmean)
        os.makedirs(racmo_arctic_data_directory+'NC_MD/'+year+'/', exist_ok=True)
        racmo_albedo_masked.to_netcdf(racmo_arctic_data_directory+'NC_MD/'+year+'/Clearsky_albedo_calculated_masked_new.nc')
        del (racmo_albedo_masked)


    if mask_albedo_constant:

        print('Masking albedo for constant date: '+constant_date)

        tilefrac5 = xr.open_dataset(racmo_arctic_data_directory + 'NC_DEFAULT/tilefrac5'+racmo_filename_additive).sel(
            time=slice(year + '-01-01', year + '-12-31'))['tilefrac5']
        tilefrac7 = xr.open_dataset(racmo_arctic_data_directory + 'NC_DEFAULT/tilefrac7'+racmo_filename_additive).sel(
            time=slice(year + '-01-01', year + '-12-31'))['tilefrac7']
        zsa = xr.open_dataset(racmo_arctic_data_directory + 'NC_DEFAULT/cosznt'+racmo_filename_additive).sel(
            time=slice(year + '-01-01', year + '-12-31'))['cosznt']

        snowcover = tilefrac5 + tilefrac7

        del(tilefrac5)
        del(tilefrac7)

        modis_albedo = xr.open_dataset(modis_data_directory + year+'_RCG.nc')['Albedo']

        modis_albedo_no = modis_albedo.where(modis_albedo['lat'] < 75)
        modis_albedo_no = modis_albedo_no.where(modis_albedo_no.values > 0.05)
        modis_albedo_no = modis_albedo_no.where(zsa.squeeze().values < np.cos(np.radians(55)))
        modis_albedo_rolmean = modis_albedo_no.rolling(time=moving_average_window, center=True).mean()

        modis_albedo_masked = modis_albedo_rolmean.where(snowcover.sel(time=year+constant_date).values == 1)
        del(modis_albedo_rolmean)
        modis_albedo_masked.to_netcdf(modis_data_directory + year + '_RCG_constant_new.nc')
        del(modis_albedo_masked)

        racmo_albedo = xr.open_dataset(racmo_arctic_data_directory+'NC_MD/Clearsky_albedo_calculated.nc').squeeze().sel(
            time=slice(year + '-01-01', year + '-12-31'))['Clear-sky_albedo']
        racmo_albedo = racmo_albedo.where(modis_albedo.values > 0.05, np.nan)
        racmo_albedo = racmo_albedo.where(racmo_albedo['lat'] <75)
        racmo_albedo = racmo_albedo.where(zsa.squeeze().values < np.cos(np.radians(55)))
        del(modis_albedo)
        racmo_albedo_rolmean = racmo_albedo.rolling(time=moving_average_window, center=True).mean()
        del(racmo_albedo)
        racmo_albedo_masked = racmo_albedo_rolmean.where(snowcover.sel(time=year+constant_date) == 1)
        del (racmo_albedo_rolmean)
        os.makedirs(racmo_arctic_data_directory+'NC_MD/'+year+'/', exist_ok=True)
        racmo_albedo_masked.squeeze().to_netcdf(racmo_arctic_data_directory+'NC_MD/'+year+'/Clearsky_albedo_calculated_constant_new.nc')
        del (racmo_albedo_masked)


    if aws_data_modification:

        t_start = year + '-01-01 00:00:00'
        t_stop = year + '-12-30 00:00:00'

        AWS = pd.read_csv(aws_directory + 'SODANKYLA_RADIATION.csv')
        AWS.rename(columns={'m': 'month', 'd': 'day'}, inplace=True)
        AWS['datetime'] = pd.to_datetime(AWS[['Year', 'month', 'day']], format="YYYY-MM-DD")
        AWS = AWS.set_index('datetime')

        AWS_daily = pd.DataFrame(index=pd.date_range(t_start, t_stop), columns=AWS.columns[4::])


        def hourly_to_daily(data):
            AWS_year_num = pd.to_numeric(data, errors='coerce')
            return (AWS_year_num.resample('D').mean())


        for i, col in enumerate(AWS_daily.columns):
            AWS_daily[col] = hourly_to_daily(AWS.loc[t_start: t_stop][col])

        AWS_daily['Reflected radiation (W/m2)'] = AWS_daily['Reflected radiation (W/m2)'] * -1
        os.makedirs(aws_directory + year + '/', exist_ok=True)
        AWS_daily.to_csv(aws_directory + year + '/SODANKYLA_RADIATION_daily.csv')

    if aws_weather:

        t_start = year + '-01-01'
        t_stop = year + '-12-30'

        AWS = pd.read_csv(aws_directory + 'SODANKYLA_WEATHER.csv')
        AWS.rename(columns={'m': 'month', 'd': 'day'}, inplace=True)
        AWS['datetime'] = pd.to_datetime(AWS[['Year', 'month', 'day']])
        AWS = AWS.set_index('datetime')

        AWS_daily = pd.DataFrame(index=np.sort(list(set(AWS.index))), columns=AWS.columns)

        for index, date in enumerate(AWS_daily.index):

            one_day_aws = AWS.loc[date,:]

            if len(one_day_aws.index) == 11:

                AWS_daily.loc[date, :] = one_day_aws
                AWS_daily.loc[date, 'Minimum temperature (degC)'] = one_day_aws.loc[
                    'Minimum temperature (degC)']
                continue

            AWS_daily.loc[one_day_aws.index[0], :] = one_day_aws.iloc[0, :]
            AWS_daily.loc[one_day_aws.index[0], 'Minimum temperature (degC)'] = one_day_aws.iloc[1,:].loc['Minimum temperature (degC)']

        AWS_daily.to_csv(aws_directory + 'SODANKYLA_WEATHER_daily.csv')

print('All years done')

if racmo_clear_sky:
    print('Calculations for clear sky radiation racmo')

    exampleds = xr.open_dataset(racmo_arctic_data_directory+'NC_MD/rsdscs'+racmo_filename_additive)
    down_shortwave = xr.open_dataset(racmo_arctic_data_directory+'NC_MD/rsdscs'+racmo_filename_additive)['rsdscs']
    net_shortwave = xr.open_dataset(racmo_arctic_data_directory+'NC_MD/ssrc'+racmo_filename_additive)['ssrc']
    up_shortwave = net_shortwave-down_shortwave

    albedo = abs(up_shortwave)/abs(down_shortwave)

    #Add metadata

    albedo = albedo.to_dataset(name='Clear-sky_albedo')

    albedo['rotated_pole'] = exampleds.rotated_pole
    albedo['Clear-sky_albedo'].attrs['Standard name'] = 'Clear sky Albedo'
    albedo['Clear-sky_albedo'].attrs['Units'] = 'W m-2'


    albedo.to_netcdf(racmo_arctic_data_directory+'NC_MD/Clearsky_albedo_calculated.nc')
    up_shortwave.to_netcdf(racmo_arctic_data_directory+'NC_MD/Shortwave_up.nc')


print('All calculations done')


