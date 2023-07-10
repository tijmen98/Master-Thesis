import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

fig_save_dir = '/Users/tijmen/Desktop/Figures_Thesis/'
statistics_dir = '/Volumes/Tijmen/Master-Thesis/Data/Statistics/'

arctic = False

norway = False

siberia = True

v1 = True
v2 = False

if v1:
    savename = 'BIAS_RMSE_snowdepth_v1.png'

if v2:
    savename = 'BIAS_RMSE_snowdepth_v2.png'

if v1 and v2:
    savename = 'BIAS_RMSE_snowdepth_v1-v2.png'


fig, axs = plt.subplots(2, 1, dpi=130, figsize=(16, 12))

years = ['2002', '2003', '2004']

month_names = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October',
               'November', 'December']

V1_BIASES = pd.DataFrame(columns=years, index = month_names)
V2_BIASES = pd.DataFrame(columns=years, index = month_names)

V1_RMSES = pd.DataFrame(columns=years, index = month_names)
V2_RMSES = pd.DataFrame(columns=years, index = month_names)

for area in ['arctic', 'norway', 'siberia']:
    arctic = False
    norway = False
    siberia = False

    if area =='arctic':
        arctic = True
        color = 'black'
    elif area =='norway':
        norway = True
        color = '#009ADE'
    elif area == 'siberia':
        siberia = True
        color = '#FF1F5B'

    for i, year in enumerate(years):

        year_dir = statistics_dir+year+'/Snowdepth/'

        if arctic:

            V1_BIAS= pd.read_csv(year_dir+'BIAS_arctic_snowdepth_v1.csv', index_col=0)
            V2_BIAS= pd.read_csv(year_dir+'BIAS_arctic_snowdepth_v2.csv', index_col=0)
            V1_RMSE= pd.read_csv(year_dir+'RMSE_arctic_snowdepth_v1.csv', index_col=0)
            V2_RMSE= pd.read_csv(year_dir+'RMSE_arctic_snowdepth_v2.csv', index_col=0)

        if norway:

            V1_BIAS= pd.read_csv(year_dir+'BIAS_norway_snowdepth_v1.csv', index_col=0)
            V2_BIAS= pd.read_csv(year_dir+'BIAS_norway_snowdepth_v2.csv', index_col=0)
            V1_RMSE= pd.read_csv(year_dir+'RMSE_norway_snowdepth_v1.csv', index_col=0)
            V2_RMSE= pd.read_csv(year_dir+'RMSE_norway_snowdepth_v2.csv', index_col=0)

        if siberia:

            V1_BIAS= pd.read_csv(year_dir+'BIAS_syberia_snowdepth_v1.csv', index_col=0)
            V2_BIAS= pd.read_csv(year_dir+'BIAS_syberia_snowdepth_v2.csv', index_col=0)
            V1_RMSE= pd.read_csv(year_dir+'RMSE_syberia_snowdepth_v1.csv', index_col=0)
            V2_RMSE= pd.read_csv(year_dir+'RMSE_syberia_snowdepth_v2.csv', index_col=0)

            V1_BIAS[6:8] = np.nan
            V2_BIAS[6:8] = np.nan
            V1_RMSE[6:8] = np.nan
            V2_RMSE[6:8] = np.nan

        V1_BIASES.loc[:,year] = list(V1_BIAS.iloc[:, 0])
        V2_BIASES.loc[:,year] = list(V2_BIAS.iloc[:, 0])
        V1_RMSES.loc[:,year] = list(V1_RMSE.iloc[:, 0])
        V2_RMSES.loc[:,year] = list(V2_RMSE.iloc[:, 0])


    test = pd.concat([pd.DataFrame(V1_BIASES.mean(axis=1)).rename(columns={0:'V1_bias'}),
                      pd.DataFrame(V2_BIASES.mean(axis=1)).rename(columns={0:'V2_bias'}),
                      pd.DataFrame(V1_RMSES.mean(axis=1)).rename(columns={0:'V1_rmse'}),
                      pd.DataFrame(V2_RMSES.mean(axis=1)).rename(columns={0:'V2_rmse'})]).to_csv(statistics_dir+savename.split('.')[0]+'.csv')

    if v1:
        axs[1].plot(np.mean(V1_RMSES.values, axis=1), color=color, linewidth=1.5, linestyle=(0, (4, 4)))
        axs[0].plot(np.mean(V1_BIASES.values, axis=1), color=color, linewidth=1.5, linestyle=(0, (4, 4)))

    if v2:
        axs[1].plot(np.mean(V2_RMSES.values, axis=1), color=color, linewidth=1.5)
        axs[0].plot(np.mean(V2_BIASES.values, axis=1), color=color, linewidth=1.5)

    if v1:
        axs[0].scatter(-0.5, np.nanmean(V1_BIASES.values), color=color, marker='v', s=20)
        axs[1].scatter(-0.5, np.nanmean(V1_RMSES.values), color=color, marker='v')

    if v2:
        axs[0].scatter(-0.5, np.nanmean(V2_BIASES.values), color=color, marker='X', s=20)
        axs[1].scatter(-0.5, np.nanmean(V2_RMSES.values), color=color, marker='X')

axs[1].set_xlim(-1, 12)
axs[1].set_xticks(np.linspace(0, 10, 6), month_names[::2])

axs[0].set_xlim(-1, 12)
axs[0].set_xticks(np.linspace(0, 10, 6), month_names[::2])

axs[0].set_ylabel('Bias [cm]')
axs[1].set_ylabel('RMSE [cm]')

axs[0].hlines(0, 0, 12, color='black', linewidth=1, zorder=-1)
axs[1].hlines(0, 0, 12, color='black', linewidth=1, zorder=-1)

legend_signs = [Line2D([0], [0], color='black', lw=4),
                Line2D([0], [0], color='#FF1F5B', lw=4),
                Line2D([0], [0], color='#009ADE', lw=4),
                Line2D([0], [0], marker='v', markersize=15, markerfacecolor='black', color='w'),
                Line2D([0], [0], marker='X', markersize=15, markerfacecolor='black', color='w'),
                Line2D([0], [0], color='black', linestyle=(0, (1, 1)), lw=4),
                Line2D([0], [0], color='black', lw=4)
                ]

if v1 and v2:
    axs[0].legend(legend_signs, ['Arctic domain',
                             'Siberia',
                             'Norway',
                             'Old scheme average value',
                             'New scheme average value',
                             'Old albedo scheme',
                             'New albedo scheme'])
elif v1:
    axs[0].legend(legend_signs[0:4], ['Arctic domain',
                             'Siberia',
                             'Norway',
                             'Average value',
                             'New scheme average value',
                             'Old albedo scheme',
                             'New albedo scheme'])

plt.savefig(fig_save_dir+savename, dpi=300)

