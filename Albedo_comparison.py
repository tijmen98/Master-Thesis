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


fig, axs = plt.subplots(3, figsize= (10,12))

year = str(2006)

WEATHER = pd.read_csv('/Volumes/Tijmen/Master-Thesis/Data/Aws_data/SODANKYLA_WEATHER_daily.csv', index_col=0).loc[
          year + '-01-01':year + '-12-31']
RADIATION = pd.read_csv('/Volumes/Tijmen/Master-Thesis/Data/Aws_data/' + year + '/SODANKYLA_RADIATION_daily.csv',
                        index_col=0)
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

    axs[1].plot(timerange, surface_albedo_to_clear_sky(THEORY_ALB_SNOWFALL(MIN_ALB, MAX_ALB, first_alb, COLD_SNOW_RELAX, timerange, snowfall_events)), color='red', linewidth=0.5)
    axs[1].plot(MODIS.loc[str(start_date):str(end_date)]['Albedo'].values, linewidth=0.5, color='orange')
    axs[1].plot(surface_albedo_to_clear_sky(AWS_ALBEDO.values), linewidth=0.5, color='green')

mean_modis_event = np.mean(modis_events, axis=0)
mean_event = np.mean(aws_events, axis=0)
mean_model_event = np.mean(model_events, axis=0)
axs[0].plot(mean_model_event, linewidth=1.5, color='blue', linestyle=(0, (1, 1)))
axs[0].plot(mean_event, linewidth=1.5, color='blue')
axs[0].plot(mean_modis_event, linewidth=1.5, color='blue', linestyle=(0, (2, 2)))

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
axs[0].plot(mean_model_event, linewidth=1.5, color='red', linestyle=(0, (1, 1)))
axs[0].plot(mean_event, linewidth=1.5, color='red')
axs[0].plot(mean_modis_event, linewidth=1.5, color='red', linestyle=(0, (2, 2)))

"""RACMO FIND MOMENTS WHERE SNOW ACCUMULATES AND EXTRACT ALBEDO TIMESERIES AROUND THESE MOMENTS"""

Precep = RACMO_SNOWFALL[RACMO_SNOWDEPTH > 2.]
racmo_events = []
model_events = []
modis_events = []
for event in range(len(Precep.index)):

    start_date = pd.to_datetime(Snowfall.index[event]) - timedelta(days=5)
    end_date = pd.to_datetime(Snowfall.index[event]) + timedelta(days=10)

    if RACMO_SNOWDEPTH[0].loc[str(start_date):str(end_date)].mean() < 10:
        continue

    snowfall = RACMO_SNOWFALL[0].loc[str(start_date):str(end_date)]
    snowfall_events = []

    for i in range(len(snowfall)):
        if snowfall.iloc[i] > 2:
            snowfall_events.append(i)

    racmo_events.append(RACMO_ALBEDO.loc[str(start_date):str(end_date)].values)
    model_events.append(THEORY_ALB_SNOWFALL(MIN_ALB, MAX_ALB, first_alb, COLD_SNOW_RELAX, timerange, snowfall_events))
    modis_events.append(MODIS.loc[str(start_date):str(end_date)]['Albedo'].values)

    axs[1].plot(timerange, surface_albedo_to_clear_sky(THEORY_ALB_SNOWFALL(MIN_ALB, MAX_ALB, first_alb, COLD_SNOW_RELAX, timerange, snowfall_events)), color='red', linewidth=0.5)
    axs[1].plot(MODIS.loc[str(start_date):str(end_date)]['Albedo'].values, linewidth=0.5, color='orange')
    axs[1].plot(surface_albedo_to_clear_sky(RACMO_ALBEDO.loc[str(start_date):str(end_date)].values), linewidth=0.5, color='green')

mean_modis_event = np.mean(modis_events, axis=0)
mean_racmo_event = np.mean(racmo_events, axis=0)
mean_model_event = np.mean(model_events, axis=0)
axs[0].plot(mean_model_event, linewidth=1.5, color='blue', linestyle=(0, (1, 1)))
axs[0].plot(mean_racmo_event, linewidth=1.5, color='blue')
axs[0].plot(mean_modis_event, linewidth=1.5, color='blue', linestyle=(0, (2, 2)))

"""RACMO FIND MOMENTS WHERE SNOW MELTS AND EXTRACT ALBEDO TIMESERIES AROUND THESE MOMENTS"""

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
axs[0].plot(mean_model_event, linewidth=1.5, color='red', linestyle=(0, (1, 1)))
axs[0].plot(mean_event, linewidth=1.5, color='red')
axs[0].plot(mean_modis_event, linewidth=1.5, color='red', linestyle=(0, (2, 2)))


axs[0].set_ylim(0,1)
axs[1].set_ylim(0,1)
axs[2].set_ylim(0,1)
plt.show()

print('Done')