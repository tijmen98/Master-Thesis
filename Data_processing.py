#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 15:02:34 2023

@author: tijmen
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime
import xarray as xr
import statistics

os.chdir('/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Python_scripts')

import Thesis_Functions.calculations as Calculations
import Thesis_Functions.data as Data


"""Variables"""

years = [2002,2003,2004]                          #list of years where data should be proccessed over, entire year is processed. Data should exist in format as specified
months = [1,2,3,4,5,6,7,8,9,10,11,12]
snow_height_threshold = 10              #Threshold in cm
breakdate = '-07-01'                    #Split date between accumulation and melt season
days_missing_limit = 5                  #Maximum number of missing days before station is discarted (MAKE MORE REFINED FILTER)

"""Calculation control"""

select_stations = True                  #Select stations that are in arctic domain
calc_stationdata = False                #Extract station snowdepth data and save to csv per station
interpolate_stationdata = False         #Interpolate
monthly_data = False                    #Extract montly data and save to directories according to structure: /Year/Month/variable.nc
monthly_data_test = True
select_stations_lat_lon = False
monthly_statistics = True

"""File names"""

racmo_snowdepth = 'NC_DEFAULT/sndp.KNMI-2001.PXARC11.RACMO24_1_complete6_UAR_q_noice_khalo6_era5q.DD.nc'


"""Directories"""

datadir = "/Volumes/Tijmen/Master-Thesis/Data/"

in_situ_data_directory = datadir+'In_situ_data/'
racmo_arctic_data_directory = datadir+'RACMO_2.4/PXARC11/2001/'
snow_cover_extend_measure_dir = datadir+'Snow_cover_Measure/'
mask_directory = datadir+'Mask/'



"""Start of yearly calculations:"""

for _ , year in enumerate(years):

    year = str(year)
    print('Starting calculations for year: '+ year)
    """Year specific directories"""
    
    in_situ_data_directory_year = in_situ_data_directory+year+'/'
    in_situ_data_directory_year_calculated = in_situ_data_directory+year+'/Calculated/'
    
    

    """If directory does not yet exists, make directories"""
    
    os.makedirs(in_situ_data_directory_year_calculated,exist_ok=True)

    """Select stations that are in arctic domain"""

    if select_stations:

        racmo_24_arc_snowheight = xr.open_dataset(racmo_arctic_data_directory + racmo_snowdepth).sel(time=slice(year + '-01-01', year + '-12-31'))['sndp']
        rlats = racmo_24_arc_snowheight['rlat'].values
        rlons = racmo_24_arc_snowheight['rlon'].values

        in_situ_data_directory_yearlist = os.listdir(in_situ_data_directory_year)
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
            observations.append(pd.read_csv(in_situ_data_directory_year + v))

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

    """Reformat raw data into station specific data, only when true"""
    
    if calc_stationdata:
        
        print('Transforming stationdata: per station yearly record')
        
        in_situ_data_directory_yearlist = os.listdir(in_situ_data_directory_year)
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
            
            observations.append(pd.read_csv(in_situ_data_directory_year+v))

        observation = pd.concat(observations)

        station_stats = pd.read_csv(in_situ_data_directory_year_calculated + 'station_in_arctic_domain_' + year + '.csv',index_col=0)

        stations = station_stats.columns

        dateindex = pd.date_range(year+"-01-01", year+"-12-31", freq="D")

        stationarray = pd.DataFrame(index=dateindex, columns =stations)

        for i,v in enumerate(stations):
            
            station_snow_depth = observation.loc[(observation['station_name'] == v) & (observation['observed_variable'] == 'snow_depth')]
            for i2 , v2 in enumerate(station_snow_depth['observation_value']):
                try: stationarray.loc[station_snow_depth['date_time'].iloc[i2].split(' ')[0],v] = v2
                
                except: stationarray.loc[station_snow_depth['date_time'].iloc[i2],v] = v2
        
        stationarray.to_csv(in_situ_data_directory_year_calculated+'stations_dayly_snowheight_'+year+'.csv')

    """Interpolate and select data with a bigger gap then days_missing_limit"""
    
    if interpolate_stationdata:
        
        print('Interpolating stationdata')
       
        stationarray = pd.read_csv(in_situ_data_directory_year_calculated+'stations_dayly_snowheight_'+year+'.csv',index_col=0)
        stations = np.array(stationarray.columns)
        
        for i,v in enumerate(stations):
    
            stationarray[v].interpolate(method='linear', inplace=True, limit = days_missing_limit)
            
            if stationarray[v].isna().sum() != 0:
                stationarray.drop(columns=[v],inplace=True)
    
        stationarray.to_csv(in_situ_data_directory_year_calculated+'stations_dayly_snowheight_interpolated_'+year+'.csv')

    """Station specific data to csv"""
    
    if monthly_data:
        
        print("saving racmo and in-situ data as monthly files, month:")
        
        """In situ data"""
        
        stationarray = pd.read_csv(in_situ_data_directory_year_calculated+'stations_dayly_snowheight_'+year+'.csv',index_col=0) 
        racmo_24_arc_snowheight = xr.open_dataset(racmo_arctic_data_directory+racmo_snowdepth).sel(time = slice(year+'-01-01',year+'-12-31'))['sndp']
        station_stats = pd.read_csv(in_situ_data_directory_year_calculated+'station_statistics_'+year+'.csv',index_col=0)
        
        startdates = pd.date_range(year+"-01-01", periods=12, freq="MS")
        enddates = pd.date_range(year+"-01-01", periods=12, freq="M")  


        
        for month in months:
            
            if len(racmo_24_arc_snowheight.time.values) < 364:
                print(racmo_24_arc_snowheight.time.values)
                print('year '+year+' is not included in the racmo dataset, please select different dataset.')
                break
            
            print(month+1)
            
            """get data for given month"""
            
            month_racmo = racmo_24_arc_snowheight.sel(time = slice(str(startdates[month]),str(enddates[month]-datetime.timedelta(days=1))))
            month_in_situ = stationarray.loc[str(startdates[month]):str(enddates[month]-datetime.timedelta(days=1))]
        
            """Check for directory to exist"""
            
            monthdir_in_situ = in_situ_data_directory_year_calculated+'/month_'+str(month+1)
            os.makedirs(monthdir_in_situ,exist_ok=True)
        
            monthdir_racmo = racmo_arctic_data_directory+year+'/month_'+str(month+1)
            os.makedirs(monthdir_racmo,exist_ok=True)
            
            """Get racmo snowheight for same locations"""
        
            locdata = pd.DataFrame(index = month_in_situ.index,columns=month_in_situ.columns)        
            
            for i,v in enumerate(station_stats.columns):
                lat = float(station_stats.loc['latitude',v])
                lon = float(station_stats.loc['longitude',v])
                
                y,x = Calculations.return_index_from_coordinates(lat,lon,racmo_24_arc_snowheight)
                
                locdata[v]=month_racmo.sel(rlat=racmo_24_arc_snowheight['rlat'][x],rlon=racmo_24_arc_snowheight['rlon'][y]).values
       
            month_in_situ.to_csv(monthdir_in_situ+'/stationdata.csv')
            locdata.to_csv(monthdir_racmo+'/stationdata.csv')
            
    if monthly_data_test:
        
        print("Test: saving racmo and in-situ data as monthly files, month:")
        
        """In situ data"""
        
        stationarray = pd.read_csv(in_situ_data_directory_year_calculated+'stations_dayly_snowheight_'+year+'.csv', index_col=0)
        racmo_24_arc_snowheight = xr.open_dataset(racmo_arctic_data_directory+racmo_snowdepth).sel(time=slice(year +'-01-01',year+'-12-31'))['sndp']
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
            
            monthdir_in_situ = in_situ_data_directory_year_calculated+'/month_'+str(month)
            os.makedirs(monthdir_in_situ,exist_ok=True)
        
            monthdir_racmo = racmo_arctic_data_directory+year+'/month_'+str(month)
            os.makedirs(monthdir_racmo,exist_ok=True)
            
            """Get racmo snowheight for same locations"""

            month_in_situ.to_csv(monthdir_in_situ+'/stationdata.csv')
            month_racmo.to_csv(monthdir_racmo+'/stationdata.csv')

    if select_stations_lat_lon == True:

        # Load the table containing the latitude and longitude information
        station_stats = pd.read_csv(in_situ_data_directory_year_calculated+'station_statistics_'+year+'.csv',index_col=0)

        
        # Define the bounding box using the minimum and maximum latitude and longitude values
        min_lat = 40.7128
        max_lat = 42.3601
        min_lon = -74.1687
        max_lon = -71.0568
        
        # Extract the locations within the bounding box
        filtered_station_stats = station_stats[(station_stats['latitude'] >= min_lat) & (station_stats['latitude'] <= max_lat) & (station_stats['longitude'] >= min_lon) & (station_stats['longitude'] <= max_lon)]
        
        filtered_station_stats.to_csv(in_situ_data_directory_year_calculated+'stations_in_latlonrange_'+year+'.csv')

    if monthly_statistics:
        station_stats = pd.read_csv(in_situ_data_directory_year_calculated+'station_in_arctic_domain_'+year+'.csv', index_col=0)

        stations = station_stats.columns
        indexes = ['racmo_mean','in_situ_mean','variance_racmo','variance_in_situ','bias']
        statistics_stations = pd.DataFrame(index=indexes,columns=stations)


        for month in months:

            print('Month ' + str(month))

            monthdir_in_situ = in_situ_data_directory_year_calculated + '/month_' + str(month)
            monthdir_racmo = racmo_arctic_data_directory + year + '/month_' + str(month)

            in_situ = pd.read_csv(monthdir_in_situ + '/stationdata.csv',index_col=0)
            racmo = pd.read_csv(monthdir_racmo + '/stationdata.csv',index_col=0)

            print(racmo.mean())

            statistics_stations['racmo_mean'] = racmo[stations].mean()
            statistics_stations['in_situ_mean'] = in_situ[stations].mean()

            statistics_stations['variance_racmo'] = racmo[stations].var()
            statistics_stations['variance_in_situ'] = in_situ[stations].var()

            statistics_stations['bias'] = racmo[stations].subtract(in_situ[stations]).mean()

            statistics_stations.to_csv(monthdir_in_situ+'/Calculated_statistics.csv')











print('All years done')