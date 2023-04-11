import os
from shapely.geometry import Polygon, Point
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib as mpl
import math
import scipy

import geopandas as gpd
import cartopy.crs as ccrs
import cartopy.feature as cfeature

os.chdir('/')

import Thesis_Functions.calculations as Calculations
import Thesis_Functions.data as Data

"""Variables"""

years = [2002, 2003, 2004]  # list of years where data should be proccessed over, entire year is processed. Data should exist in format as specified
months = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
month_names = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October',
               'November', 'December']

"""plotting control"""

arctic_domain_scatter = True
norway_scatter = True
alaska_scatter = True
canada_scatter = True
syberia_scatter = True
flat_europe_scatter = True

"""File names"""

racmo_snowdepth = 'NC_DEFAULT/sndp.KNMI-2001.PXARC11.RACMO24_1_complete6_UAR_q_noice_khalo6_era5q.DD.nc'

filename = '/Measure_merged.nc'  # Filename for combined measure dataset

"""Directories"""

fig_save_directory = '/Users/tijmen/Desktop/Figures_Thesis/'
datadir = "/Volumes/Tijmen/Master-Thesis/Data/"

in_situ_data_directory = datadir + 'In_situ_data/'
racmo_arctic_data_directory = datadir + 'RACMO_2.4/PXARC11/2001/'
snow_cover_extend_measure_dir = datadir + 'Snow_cover_Measure/'
mask_directory = datadir + 'Mask/'
remapdir = datadir + 'Remap/'
snow_cover_analysis_dir = datadir + 'Snow_cover_analyses/Snow_cover_ease_grid/'
download_measure_dir = 'Download_3-4/'

"""Bounding boxes for location extraction"""

polygon_norway = Polygon([(4.0, 55.0), (4.0, 65.0), (31.0, 75.0), (31.0, 70.0), (17.0, 66.0), (12.0, 58.0)])

polygon_flat_europe = Polygon([(4.0, 55.0), (12.0, 58.0), (17.0, 66.0), (31.0, 70.0), (31.0, 55.0)])

polygon_syberia = Polygon([(31.0, 55.0), (31.0, 75.0), (180.0, 75.0), (190.0, 66.0), (180.0, 55.0)])

polygon_alaska = Polygon([(-170.0, 55.0), (-170.0, 65.0), (-160.0, 75.0), (-160.0, 75.0), (-140.0, 55.0)])

polygon_canada = Polygon([(-130.0, 65.0), (-130.0, 85.0), (-60.0, 85.0), (-80.0, 75.0), (-60.0, 65.0)])

"""Start of yearly calculations:"""


def monthly_scatter(stations, year, racmo_directory, in_situ_directory, save_directory, save_name):

    figrange = 250

    """plot scatter heatmap day of year"""

    fig, axs = plt.subplots(3, 4, figsize=(20, 16), dpi=800)

    for month in range(12):

        xindex = 0
        yindex = month

        if yindex >= 4:
            xindex = math.floor(month / 4)
            yindex = int(month - math.floor(month / 4) * 4)

        monthdir_in_situ = in_situ_directory + 'month_' + str(month + 1)
        monthdir_racmo = racmo_directory + year + '/month_' + str(month + 1)

        in_situ = pd.read_csv(monthdir_in_situ + '/stationdata.csv', index_col=0)[stations]
        racmo = pd.read_csv(monthdir_racmo + '/stationdata.csv', index_col=0)[stations]

        in_situ = np.nan_to_num(in_situ.values, nan=-1)
        racmo = np.nan_to_num(racmo.values * 100, nan=-1)

        in_situ = in_situ.flatten('F')
        racmo = racmo.flatten('F')

        RMSE = np.sqrt(np.mean(((racmo - in_situ) ** 2)))
        try:
            regres = scipy.stats.linregress(racmo, in_situ)
        except:
            regres = np.nan
        # number_of_points = in_situ[in_situ>=0 & racmo >= 0].count()

        hist, xedges, yedges = np.histogram2d(in_situ, racmo, bins=100)
        xidx = np.clip(np.digitize(in_situ, xedges), 0, hist.shape[0] - 1)
        yidx = np.clip(np.digitize(racmo, yedges), 0, hist.shape[1] - 1)
        c = hist[xidx, yidx]

        axs[xindex, yindex].scatter(in_situ, racmo, c=c, s=1, cmap=plt.cm.RdYlBu_r, norm=mpl.colors.LogNorm())
        axs[xindex, yindex].set_title(month_names[month])
        axs[xindex, yindex].set_xlabel('In_situ')
        axs[xindex, yindex].set_ylabel('Racmo')
        axs[xindex, yindex].plot(range(figrange), range(figrange), color='black', linestyle=(0, (3, 3)), zorder=10,
                                 alpha=0.5)
        axs[xindex, yindex].set_xlim(0, figrange)
        axs[xindex, yindex].set_ylim(0, figrange)
        axs[xindex, yindex].set_aspect(1)
        axs[xindex, yindex].set_xticks(np.arange(0, figrange, figrange / 5))
        axs[xindex, yindex].set_yticks(np.arange(0, figrange, figrange / 5))
        axs[xindex, yindex].annotate(('RMSE:' + str(np.round(RMSE, 1))), xy=(150, 50))
        try:
            axs[xindex, yindex].annotate(('Slope:' + str(np.round(regres.slope, 3))), xy=(150, 40))
        except:
            axs[xindex, yindex].annotate(('Slope: none'), xy=(150, 40))
        try:
            axs[xindex, yindex].annotate(('CC:' + str(np.round(regres.rvalue, 3))), xy=(150, 30))
        except:
            axs[xindex, yindex].annotate(('Slope: none'), xy=(150, 30))

    plt.savefig(save_directory + '/' + year + '/'+save_name+'_' + year + '.png', dpi=800)


for _, year in enumerate(years):
    year = str(year)
    print('Generating plots for: ' + year)
    """Year specific directories"""

    in_situ_data_directory_year = in_situ_data_directory + year + '/'
    in_situ_data_directory_year_calculated = in_situ_data_directory + year + '/Calculated/'

    """Import area specifications"""

    station_arctic_domain = pd.read_csv(in_situ_data_directory_year_calculated + 'station_in_arctic_domain_' + year + '.csv',
                                index_col=0)

    station_stats_canada = pd.read_csv(
        in_situ_data_directory_year_calculated + 'stations_in_canada_' + year + '.csv', index_col=0)
    station_stats_syberia = pd.read_csv(
        in_situ_data_directory_year_calculated + 'stations_in_syberia_' + year + '.csv', index_col=0)
    station_stats_flat_europe = pd.read_csv(
        in_situ_data_directory_year_calculated + 'stations_in_flat_europe_' + year + '.csv', index_col=0)
    station_stats_alaska = pd.read_csv(
        in_situ_data_directory_year_calculated + 'stations_in_alaska_' + year + '.csv', index_col=0)
    station_stats_norway = pd.read_csv(
        in_situ_data_directory_year_calculated + 'stations_in_norway_' + year + '.csv', index_col=0)

    """Arctic domain scatter"""


    if arctic_domain_scatter:

        monthly_scatter(station_arctic_domain.columns.values, year, racmo_arctic_data_directory,
                        in_situ_data_directory_year_calculated, fig_save_directory, 'arctic_domain_monthly_scatter')

    if norway_scatter:

        monthly_scatter(station_stats_norway.columns.values, year, racmo_arctic_data_directory,
                        in_situ_data_directory_year_calculated, fig_save_directory, 'norway_monthly_scatter')


    if alaska_scatter:

        monthly_scatter(station_stats_alaska.columns.values, year, racmo_arctic_data_directory,
                        in_situ_data_directory_year_calculated, fig_save_directory, 'alaska_monthly_scatter')


    if canada_scatter:

        monthly_scatter(station_stats_canada.columns.values, year, racmo_arctic_data_directory,
                        in_situ_data_directory_year_calculated, fig_save_directory, 'canada_monthly_scatter')


    if flat_europe_scatter:

        monthly_scatter(station_stats_flat_europe.columns.values, year, racmo_arctic_data_directory,
                        in_situ_data_directory_year_calculated, fig_save_directory, 'flat_europe_monthly_scatter')


    if syberia_scatter:

        monthly_scatter(station_stats_syberia.columns.values, year, racmo_arctic_data_directory,
                        in_situ_data_directory_year_calculated, fig_save_directory, 'syberia_monthly_scatter')


print('all plots done')