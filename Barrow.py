import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

year = str(2003)

barrow_rad_dir = '/Volumes/Tijmen/Master-Thesis/Data/Barrow/Radiation/'
barrow_meteo_dir = '/Volumes/Tijmen/Master-Thesis/Data/Barrow/Meteo/'

cols = ['year', 'jday', 'month', 'day', 'hour', 'min', 'dt', 'zen', 'dw_solar', 'uw_solar', 'direct_n',
             'diffuse', 'dw_ir', 'dw_casetemp', 'dw_dometemp', 'uw_ir', 'uw_casetemp', 'uw_dometemp', 'uvb',
             'par', 'netsolar', 'netir', 'totalnet', 'temp', 'rh', 'windspd', 'winddir', 'pressure']


def albedo_Scheme(snowfall, temperature, old, start_alb=float(0.75), cold_relax=0.008, warm_relax=0.24, min_alb=0.5,
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


    return albedo


def read_dat_file(file_path):
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line:
                values = line.split()
                data.append(values)
    data = pd.DataFrame(data)
    return data


albedo = []

for day in range(365):

    rad_data = read_dat_file(barrow_rad_dir+year+'/brw'+format(int(year[-1]),'02d')+format(day+1, '03d')+'.dat')
    rad_data = rad_data.drop([0, 1, 2])

    dw = rad_data.iloc[:,8].astype(float)
    uw = rad_data.iloc[:,10].astype(float)

    uw = uw[dw > 1]
    dw = dw[dw > 1]

    albedo.append(np.max(uw)/np.max(dw))

albedo = np.array(albedo)

albedo[albedo > 1] = np.nan
albedo[albedo < 0] = np.nan

meteo_data = read_dat_file(barrow_meteo_dir+'met_brw_insitu_1_obop_hour_'+year+'.txt')
prec_day = meteo_data.iloc[:, 13].astype(float).groupby(meteo_data.index // 24, axis=0).sum()
temp_day = meteo_data.iloc[:, 9].astype(float).groupby(meteo_data.index // 24, axis=0).mean()

temp_day[prec_day < 0] = 0
prec_day[prec_day < 0] = 0


albedo_calc = albedo_Scheme(abs(prec_day), temp_day, old=True )
albedo_calc_new = albedo_Scheme(abs(prec_day), temp_day, old=False, cold_relax=0.15)


lims = 50, 190
xticks = np.linspace(0,365, 30)

fig, axs = plt.subplots(3)

fig.suptitle('Barrow')

axs[0].plot(albedo, color='black')
axs[0].plot(albedo_calc, color='red')
axs[0].plot(albedo_calc_new, color='blue')
axs[0].set_xticks(xticks)
axs[0].set_ylim(0, 1)
axs[0].set_xlim(lims)
axs[0].set_title('Albedo')

axs[1].plot(temp_day, color='black')
axs[1].set_xticks(xticks)
axs[1].set_ylim(-40, 20)
axs[1].set_xlim(lims)
axs[1].hlines(0, lims[0], lims[1], color='orange')
axs[1].set_title('Temperature [deg C]')

axs[2].plot(prec_day, color='black')
axs[2].set_xticks(xticks)
axs[2].set_ylim(0, 12)
axs[2].set_xlim(lims)
axs[2].set_title('Precipitation amount [mm]')

plt.tight_layout()
plt.show()

