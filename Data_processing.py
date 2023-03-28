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

os.chdir('/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Python_scripts')

import Thesis_Functions.calculations as Calculations
import Thesis_Functions.data as Data


"""Variables"""

years = [2004] #list of years where data should be proccessed over, entire year is processed. Data should exist in format as specified
months = [1,2,3,4,5,6,7,8,9,10,11]
snow_height_threshold = 10 # Threshold in cm
breakdate = '-07-01' # Split date between accumulation and melt season
days_missing_limit = 5 # Maximum number of missing days before station is discarted (MAKE MORE REFINED FILTER)

"""Calculation control"""

calc_stationdata = True            #Extract station snowdepth data and save to csv per station
interpolate_stationdata = False     #Interpolate
station_statistics = True          #Save station statistics: lat, lon
monthly_data = True                 #Extract montly data and save to directories according to structure: /Year/Month/variable.nc
monthly_data_test = False 
select_stations = False
select_stations_lat_lon = False

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
    
    """Year specific directories"""
    
    in_situ_data_directory_year = in_situ_data_directory+year+'/'
    in_situ_data_directory_year_calculated = in_situ_data_directory+year+'/Calculated/'
    
    

    """If directory does not yet exists, make directories"""
    
    os.makedirs(in_situ_data_directory_year_calculated,exist_ok=True)
    
    """reformat raw data into station specific data, only when true"""

    
    if calc_stationdata == True:
        
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
        
        stations = np.sort(list(set(observation['station_name'])))
        
        
        dateindex = pd.date_range(year+"-01-01", periods=365, freq="D")
        stationarray = pd.DataFrame(index=dateindex, columns =stations)

        for i,v in enumerate(stations):
            
            station_snow_depth = observation.loc[(observation['station_name'] == v) & (observation['observed_variable'] == 'snow_depth')]
            for i2 , v2 in enumerate(station_snow_depth['observation_value']):
                try: stationarray.loc[station_snow_depth['date_time'].iloc[i2].split(' ')[0],v] = v2
                
                except: stationarray.loc[station_snow_depth['date_time'].iloc[i2],v] = v2
        
        stationarray.to_csv(in_situ_data_directory_year_calculated+'stations_dayly_snowheight_'+year+'.csv')
    
    
    """Interpolate and select data with a bigger gap then days_missing_limit"""
    
    if interpolate_stationdata == True:
        
        print('Interpolating stationdata')
       
        stationarray = pd.read_csv(in_situ_data_directory_year_calculated+'stations_dayly_snowheight_'+year+'.csv',index_col=0)
        stations = np.array(stationarray.columns)
        
        for i,v in enumerate(stations):
    
            stationarray[v].interpolate(method='linear', inplace=True, limit = days_missing_limit)
            
            if stationarray[v].isna().sum() != 0:
                stationarray.drop(columns=[v],inplace=True)
    
        stationarray.to_csv(in_situ_data_directory_year_calculated+'stations_dayly_snowheight_interpolated_'+year+'.csv')
        
        
    """Station specific data to csv"""
    
    if station_statistics == True:
        
        print('Extracting station statistics')
        
        in_situ_data_directory_yearlist = os.listdir(in_situ_data_directory_year)
        openfiles = []
        
        for i,v in enumerate(in_situ_data_directory_yearlist):
            try: in_situ_data_directory_yearlist[i].split('.')[1] == 'csv'
            
            except:
                None
            else: 
                if in_situ_data_directory_yearlist[i].split('.')[1] == 'csv':
    
                    openfiles.append(v)
                
        observation = []
        
        for i,v in enumerate(openfiles):
            observation.append(pd.read_csv(in_situ_data_directory_year+v))
        
        observation = pd.concat(observation)
        
        stationarray = pd.read_csv(in_situ_data_directory_year_calculated+'stations_dayly_snowheight_'+year+'.csv',index_col=0)
        stations = np.array(stationarray.columns)
        indexes = ['latitude','longitude','rlat','rlon','first_day_under_threshold','last_day_under_threshold']
        
        racmo_24_arc_snowheight = xr.open_dataset(racmo_arctic_data_directory+racmo_snowdepth).sel(time = slice(year+'-01-01',year+'-12-31'))['sndp']
    
        rlats=racmo_24_arc_snowheight['rlat']
        rlons=racmo_24_arc_snowheight['rlon']
    
        station_stats = pd.DataFrame(index=indexes , columns =stations)
    
        for i,v in enumerate(stations):
            
            station_vals = observation.loc[(observation['station_name'] == v)].iloc[0]
            
            station_stats.loc['latitude',v] = station_vals['latitude']
            station_stats.loc['longitude',v] = station_vals['longitude']
            
            y,x = Calculations.return_index_from_coordinates(station_stats.loc['latitude',v],station_stats.loc['longitude',v],racmo_24_arc_snowheight)
            
            station_stats.loc['rlat',v] = rlats[x]
            station_stats.loc['rlot',v] = rlons[y]
            
            
            try: station_stats.loc['first_day_under_threshold',v] = stationarray.index[stationarray[v]>=snow_height_threshold].tolist()[0]
            except:
                station_stats.loc['first_day_under_threshold',v] = np.nan
                
            try: station_stats.loc['last_day_under_threshold',v] = stationarray.index[stationarray[v]>=snow_height_threshold].tolist()[-1]
            except:
                station_stats.loc['last_day_under_threshold',v] = np.nan
    
        station_stats.to_csv(in_situ_data_directory_year_calculated+'station_statistics_'+year+'.csv')
        

    """Extract monthly data and save to directory"""
    
    if monthly_data == True:
        
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
            
    if monthly_data_test == True:
        
        print("Test: saving racmo and in-situ data as monthly files, month:")
        
        """In situ data"""
        
        stationarray = pd.read_csv(in_situ_data_directory_year_calculated+'stations_dayly_snowheight_'+year+'.csv',index_col=0) 
        racmo_24_arc_snowheight = xr.open_dataset(racmo_arctic_data_directory+racmo_snowdepth).sel(time = slice(year+'-01-01',year+'-12-31'))['sndp']
        station_stats = pd.read_csv(in_situ_data_directory_year_calculated+'station_statistics_'+year+'.csv',index_col=0)
        
        startdates = pd.date_range(year+"-01-01", periods=12, freq="MS")
        enddates = pd.date_range(year+"-01-01", periods=12, freq="M")  
        
        rlats = station_stats.loc['rlat']
        rlons = station_stats.loc['rlon']
        
        
        for month in range(12):
            
            if len(racmo_24_arc_snowheight.time.values) != 365:
                print('year '+year+' is not included in the racmo dataset, please select different dataset.')
                break
            
            print(month+1)
            
            """get data for given month"""
            
            month_racmo = racmo_24_arc_snowheight.sel(time = slice(str(startdates[month]),str(enddates[month])),rlat=rlats,rlons=rlons)
            month_in_situ = stationarray.loc[str(startdates[month]):str(enddates[month])]
        
            """Check for directory to exist"""
            
            monthdir_in_situ = in_situ_data_directory_year_calculated+'/month_'+str(month+1)
            os.makedirs(monthdir_in_situ,exist_ok=True)
        
            monthdir_racmo = racmo_arctic_data_directory+year+'/month_'+str(month+1)
            os.makedirs(monthdir_racmo,exist_ok=True)
            
            """Get racmo snowheight for same locations"""
        
       
            month_in_situ.to_csv(monthdir_in_situ+'/stationdata.csv')
            locdata.to_csv(monthdir_racmo+'/stationdata.csv')
            
    """SELECT STATIONS"""
    
    if select_stations == True:
        
        print('Selecting stations that are in racmo range')
        
        stations_in_range = []
        
        racmo_24_arc_snowheight = xr.open_dataset(racmo_arctic_data_directory+racmo_snowdepth).sel(time = slice(year+'-01-01',year+'-12-31'))['sndp']
        station_stats = pd.read_csv(in_situ_data_directory_year_calculated+'station_statistics_'+year+'.csv',index_col=0)

        for i,v in enumerate(station_stats.columns):
            
            lat = float(station_stats.loc['latitude',v])
            lon = float(station_stats.loc['longitude',v])
            
            y,x = Calculations.return_index_from_coordinates(lat,lon,racmo_24_arc_snowheight)

            if x != 0 and x!= len(racmo_24_arc_snowheight['rlat'])-1 and y != 0 and y!= len(racmo_24_arc_snowheight['rlon'])-1:
                
                stations_in_range.append(v)
                
        pd.DataFrame(stations_in_range).to_csv(in_situ_data_directory_year_calculated+'station_in_arctic_domain_'+year+'.csv')
                
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
