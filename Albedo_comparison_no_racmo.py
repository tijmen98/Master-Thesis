"""Compare albedo time series for normal, snowfall and melting events"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import xarray as xr

"""Theoratical"""

MIN_ALB = 0.5
MAX_ALB = 0.85

COLD_SNOW_RELAX = 0.008
MELT_RELAX = 0.24

timerange = np.arange(0,15)
snowfall_event = [3, 8]
first_alb = float(0.75)


def surface_albedo_to_clear_sky(Albedo):
    # Chen et al 1983
    return(0.0587+0.730*Albedo)


def THEORY_ALB_NORMAL(MIN_ALB, MAX_ALB, start_alb, RELAX, timerange):
    albedo = np.zeros(len(timerange))
    albedo[0] = start_alb
    for time in timerange[0:-1]:
        albedo[time+1] = albedo[time] - RELAX
        if albedo[time+1] < MIN_ALB:
            albedo[time+1] = MIN_ALB
    return(albedo)
def THEORY_ALB_SNOWFALL(MIN_ALB, MAX_ALB, start_alb, RELAX, timerange, snowfall_time):
    albedo = np.zeros(len(timerange))
    albedo[0] = start_alb
    for time in timerange[0:-1]:
        albedo[time+1] = albedo[time] - RELAX
        if albedo[time+1] < MIN_ALB:
            albedo[time+1] = MIN_ALB

        if time in snowfall_time:
            albedo[time+1] = MAX_ALB
    return(albedo)
def THEORY_ALB_MELT(MIN_ALB, MAX_ALB, start_alb, RELAX, time, snowfall_time):
    albedo = np.zeros(len(timerange))
    albedo[0] = start_alb
    for time in timerange[0:-1]:
        albedo[time+1] = (albedo[time] - MIN_ALB) * np.exp(-RELAX) +MIN_ALB
        if albedo[time+1] < MIN_ALB:
            albedo[time+1] = MIN_ALB
        if time in snowfall_time:
            albedo[time+1] = MAX_ALB
    return(albedo)

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

year = str(2006)

WEATHER = pd.read_csv('/Volumes/Tijmen/Master-Thesis/Data/Sodankyla/SODANKYLA_WEATHER_daily.csv', index_col=0).loc[
          year + '-01-01':year + '-12-31']
RADIATION = pd.read_csv('/Volumes/Tijmen/Master-Thesis/Data/Sodankyla/' + year + '/SODANKYLA_RADIATION_daily.csv',
                        index_col=0)

RADIATION[abs(RADIATION) < 5] = 0

AWS_ALBEDO = abs(RADIATION['Reflected radiation (W/m2)']) / abs(
    RADIATION['Global radiation (W/m2)'])

AWS_ALBEDO[AWS_ALBEDO > 1] = np.nan
AWS_ALBEDO[AWS_ALBEDO < 0] = np.nan

lims = 200, 350
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







"""AWS AND MODIS FIND MOMENTS WHERE SNOW ACCUMULATES AND EXTRACT ALBEDO TIMESERIES AROUND THESE MOMENTS"""

Precep = WEATHER[WEATHER['Precipitation amount (mm)'] > 2.]
Snowfall = Precep[Precep['Maximum temperature (degC)'] < 0.]
aws_events = []
model_events = []
modis_events = []
for event in range(len(Snowfall.index)):

    start_date = pd.to_datetime(Snowfall.index[event]) - timedelta(days=5)
    end_date = pd.to_datetime(Snowfall.index[event]) + timedelta(days=10)

    if abs(RADIATION.loc[str(start_date):str(end_date)]['Global radiation (W/m2)']).mean() < 5:
        continue
    if WEATHER.loc[str(start_date):str(end_date)]['Snow depth (cm)'].mean() < 10:
        continue

    snowfall = WEATHER.loc[str(start_date):str(end_date)]['Precipitation amount (mm)']
    snowfall_events = []

    for i in range(len(snowfall)):
        if snowfall.iloc[i] > 2:
            snowfall_events.append(i)


    AWS_ALBEDO = abs(RADIATION.loc[str(start_date):str(end_date)]['Reflected radiation (W/m2)']) / abs(
        RADIATION.loc[str(start_date):str(end_date)]['Global radiation (W/m2)'])
    AWS_SNOWDEPTH = WEATHER.loc[str(start_date):str(end_date)]['Snow depth (cm)']

    aws_events.append(AWS_ALBEDO.values)
    model_events.append(THEORY_ALB_SNOWFALL(MIN_ALB, MAX_ALB, first_alb, COLD_SNOW_RELAX, timerange, snowfall_events))
    modis_events.append(MODIS.loc[str(start_date):str(end_date)]['Albedo'].values)

    fig, axs = plt.subplots(3, figsize=(10, 12))

    axs[1].plot(timerange, surface_albedo_to_clear_sky(THEORY_ALB_SNOWFALL(MIN_ALB, MAX_ALB, first_alb, COLD_SNOW_RELAX, timerange, snowfall_events)), color='red', linewidth=0.5)
    axs[1].plot(MODIS.loc[str(start_date):str(end_date)]['Albedo'].values, linewidth=0.5, color='orange')
    axs[1].plot(surface_albedo_to_clear_sky(AWS_ALBEDO.values), linewidth=0.5, color='green')

mean_modis_event = np.mean(modis_events, axis=0)
mean_event = np.mean(aws_events, axis=0)
mean_model_event = np.mean(model_events, axis=0)
axs[0].plot(mean_model_event, linewidth=1.5, color='blue', linestyle=(0, (1, 1)), label='mean_model_cold')
axs[0].plot(mean_event, linewidth=1.5, color='blue', label='mean_aws_cold')
axs[0].plot(mean_modis_event, linewidth=1.5, color='blue', linestyle=(0, (2, 2)), label='mean_modis_cold')

"""AWS AND MODIS FIND MOMENTS WHERE SNOW MELTS AND EXTRACT ALBEDO TIMESERIES AROUND THESE MOMENTS"""

warm = WEATHER[WEATHER['Maximum temperature (degC)'] > 0.]
aws_events = []
model_events = []
modis_events = []
for event in range(len(warm.index)):
    start_date = pd.to_datetime(warm.index[event]) - timedelta(days=5)
    end_date = pd.to_datetime(warm.index[event]) + timedelta(days=10)

    if WEATHER.loc[str(pd.to_datetime(warm.index[event])).split(' ')[0], 'Snow depth (cm)'] -2 < WEATHER.loc[str(pd.to_datetime(warm.index[event])+timedelta(days=1)).split(' ')[0], 'Snow depth (cm)']:
        continue
    if WEATHER.loc[str(start_date):str(end_date)]['Snow depth (cm)'].mean() < 10:
        continue

    snowfall = WEATHER.loc[str(start_date):str(end_date)]['Precipitation amount (mm)']

    snowfall_events = []

    for i in range(len(snowfall)):
        if snowfall.iloc[i] > 2:
            snowfall_events.append(i)

    AWS_ALBEDO = abs(RADIATION.loc[str(start_date):str(end_date)]['Reflected radiation (W/m2)']) / abs(
        RADIATION.loc[str(start_date):str(end_date)]['Global radiation (W/m2)'])
    AWS_SNOWDEPTH = WEATHER.loc[str(start_date):str(end_date)]['Snow depth (cm)']

    aws_events.append(AWS_ALBEDO.values)
    modis_events.append(MODIS.loc[str(start_date):str(end_date)]['Albedo'].values)
    model_events.append(THEORY_ALB_MELT(MIN_ALB, MAX_ALB, first_alb, MELT_RELAX, timerange, snowfall_events))

    axs[2].plot(timerange, surface_albedo_to_clear_sky(THEORY_ALB_MELT(MIN_ALB, MAX_ALB, first_alb, MELT_RELAX, timerange, snowfall_events)), color='red', linewidth=0.5)
    axs[2].plot(MODIS.loc[str(start_date):str(end_date)]['Albedo'].values, linewidth=0.5, color='orange')
    axs[2].plot(surface_albedo_to_clear_sky(AWS_ALBEDO.values), linewidth=0.5, color='green')

mean_modis_event = np.mean(modis_events, axis=0)
mean_event = np.mean(aws_events, axis=0)
mean_model_event = np.mean(model_events, axis=0)
axs[0].plot(mean_model_event, linewidth=1.5, color='red', linestyle=(0, (1, 1)), label='mean_model_warm')
axs[0].plot(mean_event, linewidth=1.5, color='red', label='mean_aws_warm')
axs[0].plot(mean_modis_event, linewidth=1.5, color='red', linestyle=(0, (2, 2)), label='mean_modis_warm')


axs[0].set_ylim(0,1)
axs[0].set_title('Mean')
axs[0].legend()
axs[1].set_ylim(0,1)
axs[1].set_title('Cold')
axs[2].set_ylim(0,1)
axs[2].set_title('Warm')

plt.show()

print('Done')