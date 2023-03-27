#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 12:11:20 2023

@author: tijmen
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

year = 2001

snow_height_threshold = 10 # Threshold in cm

breakdate = '-07-01' # Split date between accumulation and melt season

calc_stationdata = False
interpolate_stationdata = False
station_statistics = False
calculate_monthly = True

year = str(year)

datadir = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Data/In_situ_data/'+year+'/'
savedir = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Data/In_situ_data/'+year+'/Calculated/'




if calc_stationdata == True:
    
    datadirlist = os.listdir(datadir)
    openfiles = []
    
    for i,v in enumerate(datadirlist):
        if datadirlist[i].split('.')[-1] == 'csv':
            openfiles.append(v)
            
    observation = []
    
    for i,v in enumerate(openfiles):
        observation.append(pd.read_csv(datadir+v))
        
    observation = pd.concat(observation)
    
    stations = np.sort(list(set(observation['station_name'])))
    dateindex = pd.date_range("2001-01-01", periods=365, freq="D")
    stationarray = pd.DataFrame(index=dateindex , columns =stations)
    
    
    for i,v in enumerate(stations):
        
        station_snow_depth = observation.loc[(observation['station_name'] == v) & (observation['observed_variable'] == 'snow_depth')]
        for i2 , v2 in enumerate(station_snow_depth['observation_value']):
            stationarray.loc[station_snow_depth['date_time'].iloc[i2],v] = v2
    
    stationarray.to_csv(savedir+'stations_dayly_snowheight_'+year+'.csv')


if interpolate_stationdata == True:
   
    stationarray = pd.read_csv(savedir+'stations_dayly_snowheight_'+year+'.csv',index_col=0)
    stations = np.array(stationarray.columns)
    
    for i,v in enumerate(stations):

        stationarray[v].interpolate(method='linear', inplace=True, limit = 5)
        
        if stationarray[v].isna().sum() != 0:
            stationarray.drop(columns=[v],inplace=True)

    stationarray.to_csv(savedir+'stations_dayly_snowheight_interpolated_'+year+'.csv')
    
if station_statistics == True:
    
    datadirlist = os.listdir(datadir)
    openfiles = []
    
    for i,v in enumerate(datadirlist):
        if datadirlist[i].split('.')[-1] == 'csv':
            openfiles.append(v)
            
    observation = []
    
    for i,v in enumerate(openfiles):
        observation.append(pd.read_csv(datadir+v))
    
    observation = pd.concat(observation)
    
    stationarray = pd.read_csv(savedir+'stations_dayly_snowheight_interpolated_'+year+'.csv',index_col=0)
    stations = np.array(stationarray.columns)
    indexes = ['latitude','longitude','melt_season_acc_snowheight','acc_season_acc_snowheight','first_day_under_threshold','last_day_under_threshold']

    station_stats = pd.DataFrame(index=indexes , columns =stations)

    for i,v in enumerate(stations):
        
        station_vals = observation.loc[(observation['station_name'] == v)].iloc[5]
        
        station_stats.loc['latitude',v] = station_vals['latitude']
        station_stats.loc['longitude',v] = station_vals['longitude']
        
        station_stats.loc['melt_season_acc_snowheight',v] = stationarray[v].loc[year+'-01-01':year+breakdate].sum()
        station_stats.loc['acc_season_acc_snowheight',v] = stationarray[v].loc[year+breakdate:year+'-12-31'].sum()
        
        try: station_stats.loc['first_day_under_threshold',v] = stationarray.index[stationarray[v]>=snow_height_threshold].tolist()[0]
        except:
            station_stats.loc['first_day_under_threshold',v] = np.nan
            
        try: station_stats.loc['last_day_under_threshold',v] = stationarray.index[stationarray[v]>=snow_height_threshold].tolist()[-1]
        except:
            station_stats.loc['last_day_under_threshold',v] = np.nan

    station_stats.to_csv(savedir+'station_statistics_'+year+'.csv')


