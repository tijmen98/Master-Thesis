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
import xarray as xr

if laptop:
    os.chdir('/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Python_scripts')

if desktop:
    os.chdir('H:\Documenten\Master\Master_thesis\Python_scripts')
import Thesis_Functions.calculations as Calculations
import Thesis_Functions.data as Data

"""Variables"""

version = 'v2'
years = [2002]  # list of years where data should be proccessed over, entire year is processed. Data should exist in format as specified
months = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
month_names = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October',
               'November', 'December']

Albedo_masked = True
Albedo_constant = False
MODIS_filter = False

"""plotting control"""

"""Monthly scatters of snowheight in a certain domain:"""
arctic_domain_scatter = True
norway_scatter = False
alaska_scatter = False
canada_scatter = False
syberia_scatter = False
flat_europe_scatter = False

"""Map showing the study areas"""


"""Only plot for stations in tile with certain properties:"""

forest_plot = True
tundra_plot = False

area_map = False

"""File names"""

filename = '/Measure_merged.nc'  # Filename for combined measure dataset

"""Directories"""

if laptop:
    print('Directory structure: laptop')
    datadir = "/Volumes/Tijmen/Master-Thesis/Data/"
if desktop:
    print('Directory structure: desktop')
    datadir = "E:/Master-Thesis/Data/"

in_situ_data_directory = datadir + 'In_situ_data/'

if version == 'v1':
    racmo_arctic_data_directory = datadir + 'RACMO_2.4/PXARC11/2001_new/'
    racmo_filename_additive = '.KNMI-2001.PXARC11.RACMO24_1_complete6_UAR_q_noice_khalo6_era5q.DD.nc'
else:
    racmo_arctic_data_directory = datadir + 'RACMO_2.4/PXARC11/2001_v2/'
    racmo_filename_additive = '.KNMI-2001.PXARC11.RACMO24_1_complete6_UAR_q_era5q_tijmen.DD.nc'
snow_cover_extend_measure_dir = datadir + 'Snow_cover_Measure/'
mask_directory = datadir + 'Mask/'
remapdir = datadir + 'Remap/'
snow_cover_analysis_dir = datadir + 'Snow_cover_analyses/Snow_cover_ease_grid/'
download_measure_dir = 'Download_3-4/'
modis_directory = datadir+'MODIS/'
statistics_dir = '/Volumes/Tijmen/Master-Thesis/Data/Statistics/'


fig_save_directory = '/Users/tijmen/Desktop/RACMO_24_'+version+'_figures/'

"""Bounding boxes for study area"""

polygon_norway = Polygon([(4.0, 55.0), (4.0, 65.0), (31.0, 75.0), (31.0, 70.0), (17.0, 66.0), (12.0, 58.0)])

polygon_flat_europe = Polygon([(4.0, 55.0), (12.0, 58.0), (17.0, 66.0), (31.0, 70.0), (31.0, 55.0)])

polygon_syberia = Polygon([(31.0, 55.0), (31.0, 75.0), (180.0, 75.0), (190.0, 66.0), (180.0, 55.0)])

polygon_alaska = Polygon([(-170.0, 55.0), (-170.0, 65.0), (-160.0, 75.0), (-160.0, 75.0), (-140.0, 55.0)])

polygon_canada = Polygon([(-130.0, 65.0), (-130.0, 85.0), (-60.0, 85.0), (-80.0, 75.0), (-60.0, 65.0)])

"""Start of yearly calculations:"""



print('Version: '+version)

variable = 'Clear_sky_albedo'
racmo_filename = 'Clearsky_albedo_calculated_masked.nc'
savename_suffix = '_Albedo'

savename_suffix = savename_suffix+'_'+version

if tundra_plot:
    savename_suffix = savename_suffix + '_' + 'tundra'
    print('Location is Tundra')

elif forest_plot:
    print('Location is Forest')
    savename_suffix = savename_suffix + '_' + 'forest'

if Albedo_constant:
    savename_suffix = savename_suffix + '_constant_area_'


def monthly_scatter(year, var1_year, var2_year, save_directory, save_name):

    """Define limits per variable"""

    if Albedo_masked:
        limits = [0, 1]

    if Albedo_constant:
        limits = [0, 1]

    """Define labels per variable"""

    xlabel = 'Modis clear sky albedo'
    ylabel = 'Racmo clear sky albedo'


    """plot scatter heatmap day of year"""

    fig, axs = plt.subplots(3, 4, figsize=(20, 16), dpi=300)

    list_RMSE = []
    list_BIAS = []
    list_datapoints = []

    for month in range(12):

        print(month_names[month])

        xindex = 0
        yindex = month

        if yindex >= 4:
            xindex = math.floor(month / 4)
            yindex = int(month - math.floor(month / 4) * 4)

        start_date = pd.to_datetime(f'{year}-{month+1}-01')
        end_date = pd.to_datetime(f'{year}-{month+1}-{pd.Period(start_date, freq="M").days_in_month}')

        var1_month = var1_year.sel(time=slice(start_date, end_date))['Clear-sky_albedo'].squeeze()
        var2_month = var2_year.sel(time=slice(start_date, end_date))['Albedo'].squeeze()

        var1 = var1_month.values.flatten('F')
        var2 = var2_month.values.flatten('F')

        var1 = var1[var2 != 0]
        var2 = var2[var2 != 0]

        var2_nan = [not bool for bool in np.isnan(var2)]
        if var2_nan.count(True) == 0:
            continue
        var1 = var1[var2_nan]
        var2 = var2[var2_nan]

        var1_nan = [not bool for bool in np.isnan(var1)]
        if var1_nan.count(True) == 0:
            continue
        var1 = var1[var1_nan]
        var2 = var2[var1_nan]

        var2 = var2[var1 != float('inf')]
        var1 = var1[var1 != float('inf')]

        def f(B, x):
            return B[0]*x + B[1]

        linear = scipy.odr.Model(f)
        odrdata = scipy.odr.Data(var2, var1)
        odr = scipy.odr.ODR(odrdata, linear, beta0=[1, 0])
        output = odr.run()
        RMSE = np.sqrt(np.mean(((var1 - var2) ** 2)))
        bias = np.mean(var1-var2)

        list_datapoints.append(len(var1))
        list_RMSE.append(RMSE)
        list_BIAS.append(bias)

        bins = 80

        hist, xedges, yedges = np.histogram2d(var2, var1, bins=bins)
        xidx = np.clip(np.digitize(var2, xedges), 0, hist.shape[0] - 1)
        yidx = np.clip(np.digitize(var1, yedges), 0, hist.shape[1] - 1)
        c = hist[xidx, yidx]

        axs[xindex, yindex].scatter(var2, var1, c=c, s=1, cmap=plt.cm.RdYlBu_r, norm=mpl.colors.LogNorm())
        axs[xindex, yindex].set_title(month_names[month])
        axs[xindex, yindex].set_xlabel(xlabel)
        axs[xindex, yindex].set_ylabel(ylabel)
        axs[xindex, yindex].plot(np.linspace(limits[0], limits[1], 10), np.linspace(limits[0], limits[1], 10), color='black', linestyle=(0, (3, 3)), zorder=10,
                                 alpha=0.5)
        axs[xindex, yindex].plot(np.linspace(limits[0], limits[1], 10),
                                      output.beta[1] + output.beta[0]*np.linspace(limits[0], limits[1], 10),
                                      color='red', linestyle=(0, (3, 3)), zorder=20, alpha=0.5)
        axs[xindex, yindex].set_xlim(limits[0], limits[1])
        axs[xindex, yindex].set_ylim(limits[0], limits[1])
        axs[xindex, yindex].set_aspect(1)
        axs[xindex, yindex].set_xticks(np.arange(limits[0], limits[1], (limits[1]-limits[0]) / 5))
        axs[xindex, yindex].set_yticks(np.arange(limits[0], limits[1], (limits[1]-limits[0]) / 5))

        """Add statistics"""

        axs[xindex, yindex].annotate(('RMSE:' + str(np.round(RMSE, 3))), xy=(limits[0]+0.6*limits[1], limits[0]+0.25*limits[1]))
        try:
            axs[xindex, yindex].annotate(('Slope:' + str(np.round(output.beta[0], 3))), xy=(limits[0]+0.6*limits[1], limits[0]+0.3*limits[1]))
        except:
            axs[xindex, yindex].annotate(('Slope: none'), xy=(limits[0]+0.6*limits[1], limits[0]+0.3*limits[1]))
        axs[xindex, yindex].annotate(('N = '+str(len(var2))), xy=(limits[0]+0.6*limits[1], limits[0]+0.20*limits[1]))
        axs[xindex, yindex].annotate('Bias = ' + str(np.round(bias, 3)), xy=(limits[0] + 0.6 * limits[1], limits[0] + 0.15 * limits[1]))
    figure_save_directory = save_directory + '/' + year + '/' + variable + '/'
    os.makedirs(figure_save_directory, exist_ok=True)
    figure_name = save_name + '_' + year + '.png'
    plt.savefig(figure_save_directory + figure_name, dpi=300)

    os.makedirs(statistics_dir+year+'/Albedo', exist_ok=True)

    pd.DataFrame(list_RMSE).to_csv(statistics_dir+year+'/Albedo/RMSE_'+savename_suffix+'.csv')
    pd.DataFrame(list_BIAS).to_csv(statistics_dir+year+'/Albedo/BIAS_'+savename_suffix+'.csv')
    pd.DataFrame(list_datapoints).to_csv(statistics_dir + year + '/Albedo/DATAPOINTS_' + savename_suffix + '.csv')

for _, year in enumerate(years):
    year = str(year)
    print('Generating plots for: ' + year)

    """Year specific directories"""
    in_situ_data_directory_year = in_situ_data_directory + year + '/Calculated/'

    if Albedo_masked:
        var1_year = xr.open_dataset(racmo_arctic_data_directory +'NC_MD/'+ year + '/Clearsky_albedo_calculated_masked_new.nc')
        var2_year = xr.open_dataset(modis_directory+year+'_RCG_masked_new.nc')

    if Albedo_constant:
        var1_year = xr.open_dataset(racmo_arctic_data_directory +'NC_MD/'+ year + '/Clearsky_albedo_calculated_constant_new.nc')
        var2_year = xr.open_dataset(modis_directory+year+'_RCG_constant_new.nc')


    """Import area specifications"""

    if forest_plot:

        tilefrac6 = xr.open_dataset(
            racmo_arctic_data_directory + 'NC_DEFAULT/tilefrac6'+racmo_filename_additive)[
            'tilefrac6']
        tilefrac7 = xr.open_dataset(
            racmo_arctic_data_directory + 'NC_DEFAULT/tilefrac7'+racmo_filename_additive)[
            'tilefrac7']

        forest = tilefrac6.isel(time=0) + tilefrac7.isel(time=0)

        var1_year = var1_year.squeeze().where(forest.squeeze().values > 0.5)
        var2_year = var2_year.where(forest.squeeze().values > 0.5)

    if tundra_plot:

        tilefrac6 = xr.open_dataset(
            racmo_arctic_data_directory + 'NC_DEFAULT/tilefrac6'+racmo_filename_additive)[
            'tilefrac6']
        tilefrac7 = xr.open_dataset(
            racmo_arctic_data_directory + 'NC_DEFAULT/tilefrac7'+racmo_filename_additive)[
            'tilefrac7']

        forest = tilefrac6.isel(time=0) + tilefrac7.isel(time=0)

        var1_year = var1_year.squeeze().where(forest.squeeze().values < 0.1)
        var2_year = var2_year.where(forest.squeeze().values < 0.1)

    if MODIS_filter:

        var1_year = var1_year.squeeze().where(var2_year > 0.1)
        var2_year = var2_year.where(var2_year > 0.1)

    """Snowheight scatter plots"""

    if arctic_domain_scatter:

        monthly_scatter(year, var1_year,
                        var2_year, fig_save_directory, 'arctic_domain_monthly_scatter'+savename_suffix)

    if norway_scatter:

        monthly_scatter(year, var1_year,
                        var2_year, fig_save_directory, 'norway_monthly_scatter'+savename_suffix)

    if alaska_scatter:

        monthly_scatter(year, var1_year,
                        var2_year, fig_save_directory, 'alaska_monthly_scatter'+savename_suffix)

    if canada_scatter:

        monthly_scatter(year, var1_year,
                        var2_year, fig_save_directory, 'canada_monthly_scatter'+savename_suffix)

    if flat_europe_scatter:

        monthly_scatter(year, var1_year,
                        var2_year, fig_save_directory, 'flat_europe_monthly_scatter'+savename_suffix)

    if syberia_scatter:

        monthly_scatter(year, var1_year,
                        var2_year, fig_save_directory, 'syberia_monthly_scatter'+savename_suffix)

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

    plt.savefig(fig_save_directory+'map_study_areas.png', dpi=300)

print('All plots done')
