#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 14:34:01 2022

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

#%%
"""ENERGY FLUXES MAPS"""

time = 8


"""plot racmo 3 energy fluxes"""

directory3 = '/Users/tijmen/Documents/Tijmen/Climate Physics/Thesis local/Testdata/Racmo 2.3/Monthly/'
directory4 = '/Users/tijmen/Documents/Tijmen/Climate Physics/Thesis local/Testdata/Racmo 2.4/Monthly/'

directory_information = pd.read_csv(directory3+'Directory_information.csv',sep=';',names=['file','long_name','units','variable','plot'])
directory_information.drop(index=np.nan,inplace=True)
                          
sensibleheat3 = CumJoule_to_Watt(Variable_Import(directory3,'senf'))
latentheat3 = CumJoule_to_Watt(Variable_Import(directory3,'latf'))
longwaveheatflux3 = CumJoule_to_Watt(Variable_Import(directory3,'lwsn'))
shortwaveheatflux3 = CumJoule_to_Watt(Variable_Import(directory3,'swsn'))

Plot_Energy_Flux_Map(time,longwaveheatflux3,shortwaveheatflux3,sensibleheat3,latentheat3,0,mask=Get_Seamask(directory4,time))

"""plot racmo 4 energy fluxes"""
    
directory_information = pd.read_csv(directory4+'Directory_information.csv',sep=';',names=['file','long_name','units','variable','plot'])
directory_information.drop(index=np.nan,inplace=True)
                          
sensibleheat4 = CumJoule_to_Watt(Variable_Import(directory4,'hfss'))
latentheat4 = CumJoule_to_Watt(Variable_Import(directory4,'hfls'))
longwaveheatflux4 = CumJoule_to_Watt(Variable_Import(directory4,'rlds')-Variable_Import(directory4,'rlus'))
shortwaveheatflux4 = CumJoule_to_Watt(Variable_Import(directory4,'rsds')-Variable_Import(directory4,'rsus'))

Plot_Energy_Flux_Map(time,longwaveheatflux4,shortwaveheatflux4,sensibleheat4,latentheat4,0,mask=Get_Seamask(directory4,time))


"""compare RACMO 3 and EACMO 4 ENERGY fluxes"""


sensibledif = (sensibleheat4-sensibleheat3)/sensibleheat3
latentdif = (latentheat4-latentheat3)/sensibleheat3
shortwavedif = (shortwaveheatflux4-shortwaveheatflux3)/shortwaveheatflux3
longwavedif = (longwaveheatflux4-longwaveheatflux3)/longwaveheatflux3

Plot_Energy_Flux_Map_Percent(time,longwavedif,shortwavedif,sensibledif,latentdif,0,mask=Get_Seamask(directory4,time))

#%%

"""compare SMB between RACMO 3 and RACMO 4"""

directory3 = '/Users/tijmen/Documents/Tijmen/Climate Physics/Thesis local/Testdata/Racmo 2.3/Monthly/'
directory4 = '/Users/tijmen/Documents/Tijmen/Climate Physics/Thesis local/Testdata/Racmo 2.4/Monthly/'

#Racmo 3 / Racmo 4
varname = ['smb','smbgl']

time=0

Plot_Comparison_Map_Percentage(8, Variable_Import(directory3,varname[0]).isel(time=time), Variable_Import(directory4,varname[1]).isel(time=time),cmap='coolwarm',mask=Get_Seamask(directory4,time))


#%%

"""Plot number of days that tile5 is above threshold"""

directory3 = '/Users/tijmen/Documents/Tijmen/Climate Physics/Thesis local/Testdata/Racmo 2.3/Monthly/'
directory4 = '/Users/tijmen/Documents/Tijmen/Climate Physics/Thesis local/Testdata/Racmo 2.4/Monthly/'

number_of_days4 = Consecutive_Tundra(Variable_Import(directory4,'tilefrac5'),0.02)

number_of_days3 = Consecutive_Tundra(Variable_Import(directory3,'tileFR5'),0.02)

Plot_Comparison_Map_Absolute(0, number_of_days3.mean('time'), number_of_days4.mean('time'),cmap='coolwarm')


#%%



#%%

Plot_Comparison_Map_Percentage(snowmelt3.isel(time=5), snowmelt4.isel(time=5))


