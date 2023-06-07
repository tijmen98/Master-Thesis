"""Compare albedo time series for normal, snowfall and melting events"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import xarray as xr

"""Theoratical"""


def surface_albedo_to_clear_sky(Albedo):
    # Chen et al 1983
    return(0.0587+0.730*Albedo)


def albedo_Scheme(snowfall, temperature, snowdepth, old, start_alb=float(0.75), cold_relax=0.008, warm_relax=0.24, min_alb=0.5,
                  max_alb=0.85):
    timerange = range(len(snowfall.index))
    albedo = np.zeros(len(timerange))
    albedo[0] = start_alb
    for time in timerange[0:-1]:
        print(temperature[time+1])
        if temperature[time+1] < 0:
            if old:
                albedo[time+1] = albedo[time] - cold_relax
            if not old:
                albedo[time + 1] = (albedo[time] - min_alb) * np.exp(-cold_relax) + min_alb
        else:
            albedo[time + 1] = (albedo[time] - min_alb) * np.exp(-warm_relax) + min_alb

        if temperature[time+1] < 2:
            albedo[time+1] = albedo[time+1] + (np.min([np.max(snowfall[time+1], 0)/5, 1]) * (max_alb - albedo[time+1]))

        if albedo[time+1] < min_alb:
            albedo[time+1] = min_alb

        if snowdepth[time+1] < 4:
            albedo[time+1] = 0.15


    return albedo

year = str(2005)

WEATHER = pd.read_csv('/Volumes/Tijmen/Master-Thesis/Data/Sodankyla/SODANKYLA_WEATHER_daily.csv', index_col=0).loc[
          year + '-01-01':year + '-12-31']
RADIATION = pd.read_csv('/Volumes/Tijmen/Master-Thesis/Data/Sodankyla/' + year + '/SODANKYLA_RADIATION_daily.csv',
                        index_col=0)

RADIATION[abs(RADIATION) < 5] = 0

AWS_ALBEDO = abs(RADIATION['Reflected radiation (W/m2)']) / abs(
    RADIATION['Global radiation (W/m2)'])

AWS_ALBEDO[AWS_ALBEDO > 1] = np.nan
AWS_ALBEDO[AWS_ALBEDO < 0] = np.nan

lims = 50, 150
xticks = np.linspace(0,365, 12)

fig, axs = plt.subplots(3)

fig.suptitle('SODYLANKA')

axs[0].plot(albedo_Scheme(WEATHER['Precipitation amount (mm)'].astype(float), WEATHER['Air temperature (degC)'].astype(float), WEATHER['Snow depth (cm)'], old=True), color='red')
axs[0].plot(albedo_Scheme(WEATHER['Precipitation amount (mm)'].astype(float), WEATHER['Air temperature (degC)'].astype(float), WEATHER['Snow depth (cm)'], old=False, cold_relax=0.15), color='blue')
axs[0].plot(AWS_ALBEDO, color='black')
axs[0].set_xticks(xticks)
axs[0].set_ylim(0, 1)
axs[0].set_xlim(lims)
axs[0].set_title('Albedo')

axs[1].plot(WEATHER['Air temperature (degC)'], color='black')
axs[1].set_xticks(xticks)
axs[1].set_ylim(-20, 20)
axs[1].set_xlim(lims)
axs[1].hlines(0, lims[0], lims[1], color='orange')
axs[1].set_title('Temperature [deg C]')

axs[2].plot(WEATHER['Precipitation amount (mm)'], color='black')
axs[2].set_xticks(xticks)
axs[2].set_ylim(0, 12)
axs[2].set_xlim(lims)
axs[2].set_title('Precipitation amount [mm]')

plt.tight_layout()


MODIS = xr.open_dataset('/Volumes/Tijmen/Master-Thesis/Data/MODIS/'+year+'_RCG.nc').sel(rlat=26.5,
                                                                                        rlon=-11,
                                                                                        method='nearest').to_pandas()

RACMO_ALBEDO = xr.open_dataset('/Volumes/Tijmen/Master-Thesis/Data/RACMO_2.4/PXARC11/2001_new/NC_MD/Clearsky_albedo_calculated.nc')['Clear-sky_albedo'].sel(rlat=26.5,
                                                                                        rlon=-11,
                                                                                        method='nearest').to_pandas().loc[
          year + '-01-01':year + '-12-31']

RACMO_SNOWFALL = xr.open_dataset('/Volumes/Tijmen/Master-Thesis/Data/RACMO_2.4/PXARC11/2001_new/NC_DEFAULT/pr.KNMI-2001.PXARC11.RACMO24_1_complete6_UAR_q_noice_khalo6_era5q.DD.nc')['pr'].sel(rlat=26.5,
                                                                                        rlon=-11,
                                                                                        method='nearest').to_pandas().loc[
          year + '-01-01':year + '-12-31']

RACMO_SNOWDEPTH = xr.open_dataset('/Volumes/Tijmen/Master-Thesis/Data/RACMO_2.4/PXARC11/2001_new/NC_DEFAULT/sndp.KNMI-2001.PXARC11.RACMO24_1_complete6_UAR_q_noice_khalo6_era5q.DD.nc')['sndp'].sel(rlat=26.5,
                                                                                        rlon=-11,
                                                                                        method='nearest').to_pandas().loc[
          year + '-01-01':year + '-12-31']

RACMO_TEMPERATURE = xr.open_dataset('/Volumes/Tijmen/Master-Thesis/Data/RACMO_2.4/PXARC11/2001_new/NC_DEFAULT/tas.KNMI-2001.PXARC11.RACMO24_1_complete6_UAR_q_noice_khalo6_era5q.DD.nc')['tas'].sel(rlat=26.5,
                                                                                        rlon=-11,
                                                                                        method='nearest').to_pandas().loc[
          year + '-01-01':year + '-12-31']

RACMO_ALBEDO[RACMO_ALBEDO == np.inf] = np.nan