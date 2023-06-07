import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
import os


def surface_albedo_to_clear_sky(albedo):
    # Chen et al 1983
    return 0.0587+0.730*albedo


def albedo_Scheme(snowfall, temperature, start_alb=float(0.75), cold_relax=0.008, warm_relax=0.85, min_alb=0.5,
                  max_alb=0.85):
    timerange = snowfall.index
    albedo = np.zeros(len(timerange))
    albedo[0] = start_alb
    for time in timerange[0:-1]:
        if temperature[time+1] > 0:
            albedo[time+1] = albedo[time] - warm_relax
        else:
            albedo[time + 1] = (albedo[time] - min_alb) * np.exp(-cold_relax) + min_alb

        albedo[time+1] = albedo[time+1] + (np.min(np.max(snowfall[time+1], 0)/10, 1) * (max_alb - albedo[time+1]))

        if albedo[time+1] < min_alb:
            albedo[time+1] = min_alb

    return albedo

"""Define directories"""
fluxnet_main_dir = '/Volumes/Tijmen/Master-Thesis/Data/Fluxnet/'
fluxnet_data_dir = '/Volumes/Tijmen/Master-Thesis/Data/Fluxnet/Daily_fullset_Fluxnet/'

"""Open in_situ/fluxnet dir"""
ds_fluxnet_in_situ = pd.read_csv(fluxnet_main_dir+'fluxnet_in_situ.csv')

"""Import data from fluxnet"""
files = os.listdir(fluxnet_data_dir)
for file in files:
    if file[0] != 'F':
        files.remove(file)

new_files = files.copy()
print(len(new_files))

for file in files:
    file_test = file.split('_')[1]
    cont = False
    for loc in ds_fluxnet_in_situ.loc[:, 'Fluxnet'].values:
        if file.split('_')[1] == loc:
            cont = True

    if not cont:
        continue

    ds_fluxnet = pd.read_csv(fluxnet_data_dir+file, index_col=0)

    """Functions"""

    """Extract data from fluxnet dataset"""
    try:
        ds_sel = ds_fluxnet.loc[:, ['SW_IN_F_MDS', 'SW_OUT', 'P_F', 'TA_F_MDS']]
    except:
        new_files.remove(file)
        continue

    ds_sel[ds_sel == -9999] = np.nan
    ds_sel.loc['Albedo']=ds_sel['SW_OUT']/ds_sel['SW_IN_MDS']

    print(1)
