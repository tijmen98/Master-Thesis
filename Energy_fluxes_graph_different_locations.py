
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

import Thesis_Functions.calculations as calculations
import Thesis_Functions.data as data
import Thesis_Functions.plotting as plotting

"""plot graphs of energy fluxes"""

#LOAD DATA

directory3 = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/TestData/Racmo_2.3/Monthly/'
directory4 = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/TestData/Racmo_2.4/Monthly/'

directory_information = pd.read_csv(directory3+'Directory_information.csv',sep=';',names=['file','long_name','units','variable','plot'])
directory_information.drop(index=np.nan,inplace=True)
                          
sensibleheat3 = calculations.CumJoule_to_Watt_Day(data.Variable_Import(directory3,'senf'))
latentheat3 = calculations.CumJoule_to_Watt_Day(data.Variable_Import(directory3,'latf'))
longwaveheatflux3 = calculations.CumJoule_to_Watt_Day(data.Variable_Import(directory3,'lwsn'))
shortwaveheatflux3 = calculations.CumJoule_to_Watt_Day(data.Variable_Import(directory3,'swsn'))
groundheatflux3 = calculations.CumJoule_to_Watt_Day(data.Variable_Import(directory3,'gbot'))
snowmelt3 = data.Variable_Import(directory3, 'snowmelt')

directory_information = pd.read_csv(directory4+'Directory_information.csv',sep=';',names=['file','long_name','units','variable','plot'])
directory_information.drop(index=np.nan,inplace=True)
                          
sensibleheat4 = calculations.CumJoule_to_Watt_Day(data.Variable_Import(directory4,'hfss'))
latentheat4 = calculations.CumJoule_to_Watt_Day(data.Variable_Import(directory4,'hfls'))
longwaveheatflux4 = calculations.CumJoule_to_Watt_Day(data.Variable_Import(directory4,'rlds')-data.Variable_Import(directory4,'rlus'))
shortwaveheatflux4 = calculations.CumJoule_to_Watt_Day(data.Variable_Import(directory4,'rsds')-data.Variable_Import(directory4,'rsus'))
groundheatflux4 = calculations.CumJoule_to_Watt_Day(data.Variable_Import(directory4,'botflx'))
refreezing4 = calculations.CumJoule_to_Watt_Day(data.Variable_Import(directory4,'rfrzgl'))
snowmelt4 = data.Variable_Import(directory4, 'snm')

netradiation3 = calculations.Net_Radiation(shortwaveheatflux3, longwaveheatflux3, sensibleheat3, latentheat3, 0, groundheatflux3)
netradiation4 = calculations.Net_Radiation(shortwaveheatflux4, longwaveheatflux4, sensibleheat4, latentheat4, refreezing4, groundheatflux4)


areas_int = [['Canada',66.11,-72.84],
             ['Greenland',68.32,-52.1],
             ['Iceland',65.87,-22.01],
             ['Svalbard',77.99,22.77]
             ]
             


#%%

for i in range(len(areas_int)):
    fig,axs = plt.subplots(figsize=(20,16))
    

    index = calculations.return_index_from_coordinates(areas_int[i][1],areas_int[i][2],sensibleheat4)
    
    axs = plotting.Plot_Energy_Flux_Graph(index[0], index[1], longwaveheatflux4, shortwaveheatflux4, sensibleheat4, latentheat4, groundheatflux4, refreezing4,linestyle='solid',t_start = '2002-01-01',t_stop='2004-12-31')
    axs = plotting.Plot_Energy_Flux_Graph(index[0], index[1], longwaveheatflux3, shortwaveheatflux3, sensibleheat3, latentheat3,groundheatflux3, linestyle=(0,(4,4)),t_start = '2002-01-01',t_stop='2004-12-31')

    axs = plt.ylabel('Energy flux [W * m-2]')
    axs = plt.title('Comparison between RACMO3 and RACMO 4 energy fluxes in '+areas_int[i][0])
    
    
    custom_lines = [Line2D([0], [0], color='grey', lw=4,linestyle=(0,(1,1))),
                    Line2D([0], [0], color='grey', lw=4),
                Line2D([0], [0], color='Red', lw=4),
                Line2D([0], [0], color='Blue', lw=4),
                Line2D([0], [0], color='Green', lw=4),
                Line2D([0], [0], color='Yellow', lw=4),
                Line2D([0], [0], color='Orange', lw=4),
                Line2D([0], [0], color='Purple', lw=4)
                ]
    
    plt.legend(custom_lines,['RACMO 3','RACMO 4','Shortwave radiation','Longwave radiation','Sensible heatflux','Latent heatflux','Groundheatflux','Refreezing'],prop={'size': 12})
    
    plt.show()
    

#%%
    
t_start = '2002-01-01'
t_stop='2004-12-31'
      
    
for i in range(len(areas_int)):
    
    index = calculations.return_index_from_coordinates(areas_int[i][1],areas_int[i][2],sensibleheat4)
    
    fig,axs = plt.subplots(figsize=(20,16))

    snowmelt4.isel(rlat=index[1],rlon=index[0]).sel(time=slice(t_start,t_stop)).plot(label='Snowmelt RC4',color='Black',linestyle='solid',ax=axs)
    snowmelt3.isel(rlat=index[1],rlon=index[0]).sel(time=slice(t_start,t_stop)).plot(label='Snowmelt RC3',color='Black',linestyle=(0,(2,2)),ax=axs)
    
    axs2 = axs.twinx()
    
    netradiation4.isel(rlat=index[1],rlon=index[0]).sel(time=slice(t_start,t_stop)).plot(label='net radiation RC4',color='orange',linestyle='solid',ax=axs2)
    netradiation3.isel(rlat=index[1],rlon=index[0]).sel(time=slice(t_start,t_stop)).plot(label='net radiation RC3',color='orange',linestyle=(0,(2,2)),ax=axs2)
    
    axs2.set_ylabel('Heatflux [w*m-2]')
    
    custom_lines = [Line2D([0], [0], color='grey', lw=4,linestyle=(0,(1,1))),
                    Line2D([0], [0], color='grey', lw=4),
                Line2D([0], [0], color='black', lw=4),
                Line2D([0], [0], color='orange', lw=4),
                        ]
    
    plt.legend(custom_lines,['RACMO 3','RACMO 4','snowmelt','net radiation'])
        
    axs2.set_title('Comparison between RACMO3 and RACMO 4 melt and net radiation in'+areas_int[i][0])
    axs.set_title('')
    plt.show
    
    