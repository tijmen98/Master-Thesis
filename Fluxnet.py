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
import cartopy.crs as ccrs
import geopandas as gpd
import math

fluxnet_raw_to_locations = False
select_fluxnet_in_domain = False
extract_fluxnet_yearly = True

if laptop:
    os.chdir('/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Python_scripts')

if desktop:
    os.chdir('H:\Documenten\Master\Master_thesis\Python_scripts')

import Thesis_Functions.calculations as Calculations
import Thesis_Functions.data as Data

if laptop:
    print('Directory structure: laptop')
    datadir = "/Volumes/Tijmen/Master-Thesis/Data/"
if desktop:
    print('Directory structure: desktop')
    datadir = "E:/Master-Thesis/Data/"

fluxnet_data_directory = datadir+'Fluxnet/'
racmo_arctic_data_directory = datadir+'RACMO_2.4/PXARC11/2001/'


directory_list = os.listdir(fluxnet_data_directory)
fluxnet_locations = pd.read_excel(fluxnet_data_directory+'FLX_AA-Flx_BIF_DD_20200501.xlsx')



if fluxnet_raw_to_locations:

    columns = ['latitude', 'longitude']
    sites = []
    latitudes = []
    longitudes = []

    for i in fluxnet_locations.index:
        if fluxnet_locations.iloc[i, 3] == 'LOCATION_LAT':
            latitudes.append(fluxnet_locations.iloc[i, 4])
            longitudes.append(fluxnet_locations.iloc[i+1, 4])
            sites.append(fluxnet_locations.iloc[i, 0])

    data = pd.DataFrame(index=sites, columns=columns)
    data['latitude'] = latitudes
    data['longitude'] = longitudes
    data.to_csv(fluxnet_data_directory+'locations_information.csv')

if select_fluxnet_in_domain:
    locations = pd.read_csv(fluxnet_data_directory+'locations_information.csv', index_col=0)
    racmo = xr.open_dataset(racmo_arctic_data_directory+'NC_DEFAULT/tas.KNMI-2001.PXARC11.RACMO24_1_complete6_UAR_q_noice_khalo6_era5q.DD.nc')
    racmo_tas = racmo['tas'].isel(time=0).squeeze()

    rlats = racmo['rlat']
    rlons= racmo['rlon']

    stations = np.sort(locations.index)
    indexes = ['latitude', 'longitude', 'rlat', 'rlon', 'in_arctic_domain']
    station_stats = pd.DataFrame(index=indexes, columns=stations)

    print('Selecting stations that are in racmo range')

    for i, v in enumerate(stations):

        station_data = locations.loc[v]

        lat = station_data['latitude']
        lon = station_data['longitude']

        try: len(lat)
        except: None
        else:
            try: station_stats.drop(columns=v, inplace=True)
            except:
                None
            continue

        x, y = Calculations.return_index_from_coordinates(lat, lon, racmo_tas)

        rlat = rlats[y].values
        rlon = rlons[x].values

        if x != 0 and x != len(rlons) - 1 and y != 0 and y != len(rlats) - 1:
            print(v)
            station_stats.loc['in_arctic_domain', v] = True
            station_stats.loc['latitude', v] = lat
            station_stats.loc['longitude', v] = lon
            station_stats.loc['rlat', v] = rlat
            station_stats.loc['rlon', v] = rlon
        else:
            station_stats.drop(columns=v, inplace=True)

    station_stats.to_csv(fluxnet_data_directory + 'station_in_arctic_domain.csv')

if extract_fluxnet_yearly:

    stations_in_arctic_domain = pd.read_csv(fluxnet_data_directory + 'station_in_arctic_domain.csv', index_col=0)

    for i, v in enumerate(directory_list):
        if v[-1] != 'p' and v[0] != '.':

            string = v.split('_')
            try: string[4]
            except: continue
            else:
                if string[4] == 'DD' and string[1] in stations_in_arctic_domain.columns:
                    locationdata = pd.read_csv(fluxnet_data_directory+v)
                    print(locationdata)

