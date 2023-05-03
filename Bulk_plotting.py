desktop = False
laptop = True

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

if laptop:
    os.chdir('/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Python_scripts')

if desktop:
    os.chdir('H:\Documenten\Master\Master_thesis\Python_scripts')
import Thesis_Functions.calculations as Calculations
import Thesis_Functions.data as Data

"""Variables"""

years = [2002, 2003, 2004]  # list of years where data should be proccessed over, entire year is processed. Data should exist in format as specified
months = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
month_names = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October',
               'November', 'December']

Snowdepth = True
Surface_temp = False
Precipitation = False

"""plotting control"""

"""Data type"""

no_nan_data = False

"""Monthly scatters of snowheight in a certain domain:"""
arctic_domain_scatter = True
norway_scatter = False
alaska_scatter = False
canada_scatter = False
syberia_scatter = False
flat_europe_scatter = False

"""Map showing the study areas"""


"""Only map stations in tile with certain tilefraction"""

tilefractionplotting = True
tilefrac = 'tilefrac9'

area_map = False

"""File names"""

racmo_snowdepth = 'NC_DEFAULT/sndp.KNMI-2001.PXARC11.RACMO24_1_complete6_UAR_q_noice_khalo6_era5q.DD.nc'
filename = '/Measure_merged.nc'  # Filename for combined measure dataset

"""Directories"""

fig_save_directory = '/Users/tijmen/Desktop/Figures_Thesis/'

if laptop:
    print('Directory structure: laptop')
    datadir = "/Volumes/Tijmen/Master-Thesis/Data/"
if desktop:
    print('Directory structure: desktop')
    datadir = "E:/Master-Thesis/Data/"

in_situ_data_directory = datadir + 'In_situ_data/'
racmo_arctic_data_directory = datadir + 'RACMO_2.4/PXARC11/2001/'
snow_cover_extend_measure_dir = datadir + 'Snow_cover_Measure/'
mask_directory = datadir + 'Mask/'
remapdir = datadir + 'Remap/'
snow_cover_analysis_dir = datadir + 'Snow_cover_analyses/Snow_cover_ease_grid/'
download_measure_dir = 'Download_3-4/'

"""Bounding boxes for study area"""

polygon_norway = Polygon([(4.0, 55.0), (4.0, 65.0), (31.0, 75.0), (31.0, 70.0), (17.0, 66.0), (12.0, 58.0)])

polygon_flat_europe = Polygon([(4.0, 55.0), (12.0, 58.0), (17.0, 66.0), (31.0, 70.0), (31.0, 55.0)])

polygon_syberia = Polygon([(31.0, 55.0), (31.0, 75.0), (180.0, 75.0), (190.0, 66.0), (180.0, 55.0)])

polygon_alaska = Polygon([(-170.0, 55.0), (-170.0, 65.0), (-160.0, 75.0), (-160.0, 75.0), (-140.0, 55.0)])

polygon_canada = Polygon([(-130.0, 65.0), (-130.0, 85.0), (-60.0, 85.0), (-80.0, 75.0), (-60.0, 65.0)])

"""Start of yearly calculations:"""

if Snowdepth:
    print('Variable: Snowdepth')

    in_situ_variable = 'snow_depth'
    racmo_filename = 'NC_DEFAULT/sndp.KNMI-2001.PXARC11.RACMO24_1_complete6_UAR_q_noice_khalo6_era5q.DD.nc'
    racmo_variable = 'sndp'
    savename_suffix = 'snowdepth'

if Surface_temp:
    print('Variable: Surface temperature')

    in_situ_variable = 'air_temperature'
    racmo_filename = 'NC_DEFAULT/tas.KNMI-2001.PXARC11.RACMO24_1_complete6_UAR_q_noice_khalo6_era5q.DD.nc'
    racmo_variable = 'tas'
    savename_suffix = 'surface_temperature'

if Precipitation:
    print('Variable: Precipitation')

    in_situ_variable = 'accumulated_precipitation'
    racmo_filename = 'NC_DEFAULT/pr.KNMI-2001.PXARC11.RACMO24_1_complete6_UAR_q_noice_khalo6_era5q.DD.nc'
    racmo_variable = 'pr'
    savename_suffix = 'precipitation'

if tilefractionplotting:
    savename_suffix = savename_suffix + '_' + tilefrac
def monthly_scatter(stations, year, racmo_directory, in_situ_directory, save_directory, save_name):

    if Snowdepth:
        limits = [-1, 150]
    if Surface_temp:
        limits = [230, 330]

    """plot scatter heatmap day of year"""

    fig, axs = plt.subplots(3, 4, figsize=(20, 16), dpi=800)

    for month in range(12):

        xindex = 0
        yindex = month

        if yindex >= 4:
            xindex = math.floor(month / 4)
            yindex = int(month - math.floor(month / 4) * 4)

        monthdir_in_situ = in_situ_directory + 'month_' + str(month + 1)
        monthdir_racmo = racmo_directory + '/month_' + str(month + 1)

        in_situ = pd.read_csv(monthdir_in_situ + '/stationdata.csv', index_col=0)[stations]
        racmo = pd.read_csv(monthdir_racmo + '/stationdata.csv', index_col=0)[stations]

        in_situ = in_situ.values.flatten('F')
        racmo = racmo.values.flatten('F')

        in_situ_nan = [not bool for bool in np.isnan(in_situ)]
        if not Snowdepth:
            racmo = racmo[in_situ_nan]
        if Snowdepth:
            racmo = racmo[in_situ_nan]*100
        in_situ = in_situ[in_situ_nan]

        def f(B, x):
            return B[0]*x + B[1]

        linear = scipy.odr.Model(f)
        odrdata = scipy.odr.Data(in_situ, racmo)
        odr = scipy.odr.ODR(odrdata, linear, beta0=[0.5, 0])
        output = odr.run()

        RMSE = np.sqrt(np.mean(((racmo - in_situ) ** 2)))

        bins = np.round((len(in_situ)/1000)*4,0).astype(int)
        if bins < 15:
            bins = int(15)

        hist, xedges, yedges = np.histogram2d(in_situ, racmo, bins=bins)
        xidx = np.clip(np.digitize(in_situ, xedges), 0, hist.shape[0] - 1)
        yidx = np.clip(np.digitize(racmo, yedges), 0, hist.shape[1] - 1)
        c = hist[xidx, yidx]

        axs[xindex, yindex].scatter(in_situ, racmo, c=c, s=1, cmap=plt.cm.RdYlBu_r, norm=mpl.colors.LogNorm())
        axs[xindex, yindex].set_title(month_names[month])
        axs[xindex, yindex].set_xlabel('In-situ snowheight [cm]')
        axs[xindex, yindex].set_ylabel('Racmo snowheight [cm]')
        axs[xindex, yindex].plot(np.arange(limits[0], limits[1]), np.arange(limits[0], limits[1]), color='black', linestyle=(0, (3, 3)), zorder=10,
                                 alpha=0.5)
        try: axs[xindex, yindex].plot(np.arange(limits[0], limits[1]), output.beta[1] + output.beta[0]*np.arange(limits[0], limits[1]), color='red',
                                      linestyle=(0, (3, 3)), zorder=10, alpha=0.5)
        except: None
        axs[xindex, yindex].set_xlim(limits[0], limits[1])
        axs[xindex, yindex].set_ylim(limits[0], limits[1])
        axs[xindex, yindex].set_aspect(1)
        axs[xindex, yindex].set_xticks(np.arange(limits[0], limits[1], (limits[1]-limits[0]) / 5))
        axs[xindex, yindex].set_yticks(np.arange(limits[0], limits[1], (limits[1]-limits[0]) / 5))
        axs[xindex, yindex].annotate(('RMSE:' + str(np.round(RMSE, 1))), xy=(limits[1]-30, limits[0]+5*((limits[1]-limits[0])/20)))
        try:
            axs[xindex, yindex].annotate(('Slope:' + str(np.round(output.beta[0], 3))), xy=(limits[1]-30, limits[0]+4*((limits[1]-limits[0])/20)))
        except:
            axs[xindex, yindex].annotate(('Slope: none'), xy=(limits[1]-30, limits[0]+4*((limits[1]-limits[0])/20)))
        try:
            axs[xindex, yindex].annotate(('CC:' + str(np.round(regres.rvalue, 3))), xy=(limits[1]-30, limits[0]+3*((limits[1]-limits[0])/20)))
        except:
            axs[xindex, yindex].annotate(('CC: none'), xy=(limits[1]-30, limits[0]+3*((limits[1]-limits[0])/20)))
        axs[xindex, yindex].annotate(('N = '+str(len(in_situ))), xy=(limits[1]-30, limits[0]+2*((limits[1]-limits[0])/20)))

    plt.savefig(save_directory + '/' + year + '/'+save_name+'_' + year + '.png', dpi=800)

for _, year in enumerate(years):
    year = str(year)
    print('Generating plots for: ' + year)
    """Year specific directories"""

    in_situ_data_directory_year = in_situ_data_directory + year + '/Calculated/'
    in_situ_data_directory_year_calculated = in_situ_data_directory + year + '/Calculated/'+in_situ_variable+'/'

    if no_nan_data:
        in_situ_data_directory_year_calculated = in_situ_data_directory + year + '/Calculated/' + in_situ_variable + '_no_nan/'

    racmo_arctic_data_directory_year = racmo_arctic_data_directory + racmo_variable +'/'+ year

    """Import area specifications"""

    station_arctic_domain = pd.read_csv(in_situ_data_directory_year + '/station_in_arctic_domain_' + year + '.csv',
                                index_col=0)

    station_stats_canada = pd.read_csv(
        in_situ_data_directory_year + 'stations_in_canada_' + year + '.csv', index_col=0)
    station_stats_syberia = pd.read_csv(
        in_situ_data_directory_year + 'stations_in_syberia_' + year + '.csv', index_col=0)
    station_stats_flat_europe = pd.read_csv(
        in_situ_data_directory_year + 'stations_in_flat_europe_' + year + '.csv', index_col=0)
    station_stats_alaska = pd.read_csv(
        in_situ_data_directory_year + 'stations_in_alaska_' + year + '.csv', index_col=0)
    station_stats_norway = pd.read_csv(
        in_situ_data_directory_year + 'stations_in_norway_' + year + '.csv', index_col=0)

    if tilefractionplotting:
        station_arctic_domain = pd.read_csv(
            in_situ_data_directory_year + 'station_'+tilefrac + '_' + year + '.csv', index_col=0)

    """Snowheight scatter plots"""

    if arctic_domain_scatter:

        monthly_scatter(station_arctic_domain.columns.values, year, racmo_arctic_data_directory_year,
                        in_situ_data_directory_year_calculated, fig_save_directory, 'arctic_domain_monthly_scatter'+savename_suffix)

    if norway_scatter:

        monthly_scatter(station_stats_norway.columns.values, year, racmo_arctic_data_directory_year,
                        in_situ_data_directory_year_calculated, fig_save_directory, 'norway_monthly_scatter'+savename_suffix)

    if alaska_scatter:

        monthly_scatter(station_stats_alaska.columns.values, year, racmo_arctic_data_directory_year,
                        in_situ_data_directory_year_calculated, fig_save_directory, 'alaska_monthly_scatter'+savename_suffix)

    if canada_scatter:

        monthly_scatter(station_stats_canada.columns.values, year, racmo_arctic_data_directory_year,
                        in_situ_data_directory_year_calculated, fig_save_directory, 'canada_monthly_scatter'+savename_suffix)

    if flat_europe_scatter:

        monthly_scatter(station_stats_flat_europe.columns.values, year, racmo_arctic_data_directory_year,
                        in_situ_data_directory_year_calculated, fig_save_directory, 'flat_europe_monthly_scatter'+savename_suffix)

    if syberia_scatter:

        monthly_scatter(station_stats_syberia.columns.values, year, racmo_arctic_data_directory_year,
                        in_situ_data_directory_year_calculated, fig_save_directory, 'syberia_monthly_scatter'+savename_suffix)

    """Snow extend scatter plots"""

"""Non yearly plots:"""

print('Generating non yearly plots:')

if area_map:

    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={'projection': ccrs.Stereographic(central_longitude=0.,
                                                                                        central_latitude=90.)})
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS)
    gpd.GeoDataFrame(geometry=[polygon_norway]).plot(ax=ax, transform=ccrs.PlateCarree(), color='red')
    gpd.GeoDataFrame(geometry=[polygon_flat_europe]).plot(ax=ax, transform=ccrs.PlateCarree(), color='green')
    gpd.GeoDataFrame(geometry=[polygon_syberia]).plot(ax=ax, transform=ccrs.PlateCarree(), color='blue')
    gpd.GeoDataFrame(geometry=[polygon_alaska]).plot(ax=ax, transform=ccrs.PlateCarree(), color='orange')
    gpd.GeoDataFrame(geometry=[polygon_canada]).plot(ax=ax, transform=ccrs.PlateCarree(), color='yellow')

    plt.savefig(fig_save_directory+'map_study_areas.png', dpi=500)

print('All plots done')
