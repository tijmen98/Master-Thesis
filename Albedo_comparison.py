"""Compare albedo time series for normal, snowfall and melting events"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

"""Theoratical"""

MIN_ALB = 0.5
MAX_ALB = 0.85

COLD_SNOW_RELAX = 0.008
MELT_RELAX = 0.24

timerange = np.arange(0,60)
snowfall_event = 25
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
def THEORY_ALB_MELT(MIN_ALB, MAX_ALB, start_alb, RELAX, time):
    albedo = np.zeros(len(timerange))
    albedo[0] = start_alb
    for time in timerange[0:-1]:
        albedo[time+1] = (albedo[time] - MIN_ALB) * np.exp(-RELAX) +MIN_ALB
        if albedo[time+1] < MIN_ALB:
            albedo[time+1] = MIN_ALB
    return(albedo)

plt.plot(timerange,THEORY_ALB_NORMAL(MIN_ALB, MAX_ALB, first_alb, COLD_SNOW_RELAX, timerange))
plt.plot(timerange,THEORY_ALB_SNOWFALL(MIN_ALB, MAX_ALB, first_alb, COLD_SNOW_RELAX, timerange, snowfall_event))
plt.plot(timerange,THEORY_ALB_MELT(MIN_ALB, MAX_ALB, first_alb, MELT_RELAX, timerange))
plt.show()

print('Done')