"""Compare albedo time series for normal, snowfall and melting events"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

"""Theoratical"""

MIN_ALB = 0.5
MAX_ALB = 0.85

COLD_SNOW_RELAX = 0.008
MELT_RELAX = 0.24

timerange = np.arange(0,15)
snowfall_event = 3
first_alb = float(0.75)


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

        if time == snowfall_time:
            albedo[time+1] = MAX_ALB
    return(albedo)
def THEORY_ALB_MELT(MIN_ALB, MAX_ALB, start_alb, RELAX, time, snowfall_time):
    albedo = np.zeros(len(timerange))
    albedo[0] = start_alb
    for time in timerange[0:-1]:
        albedo[time+1] = (albedo[time] - MIN_ALB) * np.exp(-RELAX) +MIN_ALB
        if albedo[time+1] < MIN_ALB:
            albedo[time+1] = MIN_ALB
        if time == snowfall_time:
            albedo[time+1] = MAX_ALB
    return(albedo)


fig, axs = plt.subplots(3, figsize= (10,12))

axs[0].plot(timerange,THEORY_ALB_SNOWFALL(MIN_ALB, MAX_ALB, first_alb, COLD_SNOW_RELAX, timerange, snowfall_event))
axs[0].plot(timerange,THEORY_ALB_MELT(MIN_ALB, MAX_ALB, first_alb, MELT_RELAX, timerange, snowfall_event))

year= str(2005)

WEATHER = pd.read_csv('/Volumes/Tijmen/Master-Thesis/Data/Aws_data/SODANKYLA_WEATHER_daily.csv', index_col=0).loc[
          year + '-01-01':year + '-12-31']
RADIATION = pd.read_csv('/Volumes/Tijmen/Master-Thesis/Data/Aws_data/' + year + '/SODANKYLA_RADIATION_daily.csv',
                        index_col=0)

"""FIND MOMENTS WHERE SNOW ACCUMULATES AND EXTRACT ALBEDO TIMESERIES AROUND THESE MOMENTS"""

Precep = WEATHER[WEATHER['Precipitation amount (mm)'] > 2.]
Snowfall = Precep[Precep['Maximum temperature (degC)'] < 0.]

events = []

for event in range(len(Snowfall.index)):
    start_date = pd.to_datetime(Snowfall.index[event]) - timedelta(days=5)
    end_date = pd.to_datetime(Snowfall.index[event]) + timedelta(days=10)

    if abs(RADIATION.loc[str(start_date):str(end_date)]['Global radiation (W/m2)']).mean() < 5:
        continue
    if WEATHER.loc[str(start_date):str(end_date)]['Snow depth (cm)'].mean() < 10:
        continue

    AWS_ALBEDO = abs(RADIATION.loc[str(start_date):str(end_date)]['Reflected radiation (W/m2)']) / abs(
        RADIATION.loc[str(start_date):str(end_date)]['Global radiation (W/m2)'])
    AWS_SNOWDEPTH = WEATHER.loc[str(start_date):str(end_date)]['Snow depth (cm)']
    events.append(AWS_ALBEDO.values)
    axs[1].plot(AWS_ALBEDO.values)
mean_event = np.mean(events, axis=0)

axs[1].plot(mean_event, linewidth=3)
axs[0].plot(mean_event, linewidth=3)

axs[0].set_ylim(0,1)
axs[1].set_ylim(0,1)
plt.show()

"""FIND MOMENTS WERE SNOW MELTS AND EXTRACT ALBEDO TIMESERIES AROUND THESE MOMENTS"""

print('Done')