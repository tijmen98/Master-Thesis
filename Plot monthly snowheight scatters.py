#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 14:45:15 2023

@author: tijmen
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib as mpl
import math
import scipy

year = str(2004)


datadir = "/Volumes/Tijmen/Master-Thesis/Data/"

in_situ_data_directory = datadir+'In_situ_data/'
in_situ_data_directory_year_calculated = in_situ_data_directory+year+'/Calculated/'
racmo_arctic_data_directory = datadir+'RACMO_2.4/PXARC11/2001/'

fig_save_directory = '/Users/tijmen/Desktop/Figures_Thesis/'

scatterplot = True

figrange = 250

stationselect = pd.read_csv(in_situ_data_directory_year_calculated+'station_in_arctic_domain_'+year+'.csv',index_col=0)

if scatterplot == True:
    
    """plot scatter heatmap day of year"""
    
    fig, axs = plt.subplots(3,4,figsize=(20,16),dpi=800)
    
    month_names = ['January','February','March','April','May','June','July','August','September','October','November','December']
    
    for month in range(12):
        
        xindex = 0
        yindex = month
        
        if yindex >= 4:
            xindex = math.floor(month/4)
            yindex = int(month - math.floor(month/4)*4)

        monthdir_in_situ = in_situ_data_directory_year_calculated+'month_'+str(month+1)
        monthdir_racmo = racmo_arctic_data_directory+year+'/month_'+str(month+1)

        in_situ = pd.read_csv(monthdir_in_situ+'/stationdata.csv',index_col=0).loc[:,stationselect.columns]
        racmo = pd.read_csv(monthdir_racmo+'/stationdata.csv',index_col=0).loc[:,stationselect.columns]
        
        in_situ = np.nan_to_num(in_situ.values,nan=-1)
        racmo = np.nan_to_num(racmo.values*100,nan=-1)
        
        in_situ = in_situ.flatten('F')
        racmo = racmo.flatten('F')
        
        RMSE = np.sqrt(np.mean(((racmo-in_situ)**2)))
        regres =scipy.stats.linregress(racmo,in_situ)
        #number_of_points = in_situ[in_situ>=0 & racmo >= 0].count()
    
        hist, xedges, yedges = np.histogram2d(in_situ,racmo,bins=100)
        xidx = np.clip(np.digitize(in_situ, xedges), 0, hist.shape[0]-1)
        yidx = np.clip(np.digitize(racmo, yedges), 0, hist.shape[1]-1)
        c = hist[xidx, yidx]
        
    
        axs[xindex,yindex].scatter(in_situ, racmo, c=c, s = 1, cmap=plt.cm.RdYlBu_r,norm=mpl.colors.LogNorm())
        axs[xindex,yindex].set_title(month_names[month])
        axs[xindex,yindex].set_xlabel('In_situ')
        axs[xindex,yindex].set_ylabel('Racmo')
        axs[xindex,yindex].plot(range(figrange),range(figrange),color='black',linestyle=(0,(3,3)),zorder=10,alpha=0.5)
        axs[xindex,yindex].set_xlim(0,figrange)
        axs[xindex,yindex].set_ylim(0,figrange)
        axs[xindex,yindex].set_aspect(1)
        axs[xindex,yindex].set_xticks(np.arange(0,figrange,figrange/5))
        axs[xindex,yindex].set_yticks(np.arange(0,figrange,figrange/5))
        axs[xindex,yindex].annotate(('RMSE:'+str(np.round(RMSE,1))),xy=(150,50))
        axs[xindex,yindex].annotate(('Slope:'+str(np.round(regres.slope,3))),xy=(150,40))  
        axs[xindex,yindex].annotate(('CC:'+str(np.round(regres.rvalue,3))),xy=(150,30))
    
    
    plt.savefig(fig_save_directory+'/'+year+'/RACMO_monthly_scatter_'+year+'.png',dpi=800)
    