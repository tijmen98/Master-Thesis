import numpy as np
import matplotlib.pyplot as plt

fig_save_dir = '/Users/tijmen/Desktop/Figures_Thesis/'

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

fig_save_dir = '/Users/tijmen/Desktop/Figures_Thesis/'
statistics_dir = '/Volumes/Tijmen/Master-Thesis/Data/Statistics/'

v1 = True
v2 = False

if v1:
    savename = 'BIAS_RMSE_albedo_v1.png'

if v2:
    savename = 'BIAS_RMSE_albedo_v2.png'

if v1 and v2:
    savename = 'BIAS_RMSE_albedo_v1-v2.png'


fig, axs = plt.subplots(3, 1, dpi=130, figsize=(16, 12), gridspec_kw={'height_ratios': [2, 2, 1]}, sharex=True)
axs[2].set_yscale('linear')
axs[2].annotate('Maximum number of datapoints:', xy=(0, 27), size=5)

for x in range(3):
    axs[x].grid()

years = ['2002', '2003', '2004', '2005']

month_names = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October',
               'November', 'December']

V1_BIASES = pd.DataFrame(columns=years, index = month_names)
V2_BIASES = pd.DataFrame(columns=years, index = month_names)

V1_RMSES = pd.DataFrame(columns=years, index = month_names)
V2_RMSES = pd.DataFrame(columns=years, index = month_names)

V1_DPS = pd.DataFrame(columns=years, index = month_names)
V2_DPS = pd.DataFrame(columns=years, index = month_names)

for i_area, area in enumerate(['Arctic', 'Constant', 'Tundra', 'Forest']):
    arctic = False
    constant = False
    tundra = False
    forest = False

    print(area)

    if area == 'Arctic':
        arctic = True
        color = 'black'
    elif area == 'Constant':
        constant = True
        color = '#A0B1BA'
    elif area == 'Forest':
        forest = True
        color = '#00CD6C'
    elif area == 'Tundra':
        tundra = True
        color = '#A6761D'
    
    for i, year in enumerate(years):
    
        year_dir = statistics_dir+year+'/Albedo/'
    
        if arctic:

    
            V1_BIAS= pd.read_csv(year_dir+'BIAS__Albedo_v1.csv', index_col=0)
            V2_BIAS= pd.read_csv(year_dir+'BIAS__Albedo_v2.csv', index_col=0)
            V1_RMSE= pd.read_csv(year_dir+'RMSE__Albedo_v1.csv', index_col=0)
            V2_RMSE= pd.read_csv(year_dir+'RMSE__Albedo_v2.csv', index_col=0)
            V1_DP = pd.read_csv(year_dir+'DATAPOINTS__Albedo_v1.csv', index_col=0)
            V2_DP = pd.read_csv(year_dir+'DATAPOINTS__Albedo_v1.csv', index_col=0)
    
        if constant:
    
            V1_BIAS= pd.read_csv(year_dir+'BIAS__Albedo_v1_constant_area_.csv', index_col=0)
            V2_BIAS= pd.read_csv(year_dir+'BIAS__Albedo_v2_constant_area_.csv', index_col=0)
            V1_RMSE= pd.read_csv(year_dir+'RMSE__Albedo_v1_constant_area_.csv', index_col=0)
            V2_RMSE= pd.read_csv(year_dir+'RMSE__Albedo_v2_constant_area_.csv', index_col=0)
            V1_DP = pd.read_csv(year_dir+'DATAPOINTS__Albedo_v1_constant_area_.csv', index_col=0)
            V2_DP = pd.read_csv(year_dir+'DATAPOINTS__Albedo_v2_constant_area_.csv', index_col=0)

        if tundra:

            V1_BIAS= pd.read_csv(year_dir+'BIAS__Albedo_v1_tundra.csv', index_col=0)
            V2_BIAS= pd.read_csv(year_dir+'BIAS__Albedo_v2_tundra.csv', index_col=0)
            V1_RMSE= pd.read_csv(year_dir+'RMSE__Albedo_v1_tundra.csv', index_col=0)
            V2_RMSE= pd.read_csv(year_dir+'RMSE__Albedo_v2_tundra.csv', index_col=0)
            V1_DP = pd.read_csv(year_dir+'DATAPOINTS__Albedo_v1_tundra.csv', index_col=0)
            V2_DP = pd.read_csv(year_dir+'DATAPOINTS__Albedo_v1_tundra.csv', index_col=0)


        if forest:

            V1_BIAS= pd.read_csv(year_dir+'BIAS__Albedo_v1_forest.csv', index_col=0)
            V2_BIAS= pd.read_csv(year_dir+'BIAS__Albedo_v2_forest.csv', index_col=0)
            V1_RMSE= pd.read_csv(year_dir+'RMSE__Albedo_v1_forest.csv', index_col=0)
            V2_RMSE= pd.read_csv(year_dir+'RMSE__Albedo_v2_forest.csv', index_col=0)
            V1_DP = pd.read_csv(year_dir+'DATAPOINTS__Albedo_v1_forest.csv', index_col=0)
            V2_DP = pd.read_csv(year_dir+'DATAPOINTS__Albedo_v1_forest.csv', index_col=0)
    
    
        V1_BIASES.loc[:,year] = list(V1_BIAS.iloc[:, 0])
        V2_BIASES.loc[:,year] = list(V2_BIAS.iloc[:, 0])
        V1_RMSES.loc[:,year] = list(V1_RMSE.iloc[:, 0])
        V2_RMSES.loc[:,year] = list(V2_RMSE.iloc[:, 0])
        V1_DPS.loc[:,year] = list(V1_DP.iloc[:, 0])
        V2_DPS.loc[:,year] = list(V2_DP.iloc[:, 0])

    axs[2].plot(np.nanmean((V1_DPS/np.max(V1_DPS)*100), axis=1), color=color, linewidth=1)
    axs[2].annotate(area+': '+str(int(np.round(np.mean(np.max(V1_DPS)), 0))), xy=(0, 20 - i_area*7), size=5)

    
    pd.concat([pd.DataFrame(V1_BIASES),
                      pd.DataFrame(V2_BIASES)], keys=['old', 'new']).to_csv(statistics_dir+savename.split('.')[0]+'_BIAS.csv', decimal=',')
    
    pd.concat([pd.DataFrame(V1_RMSES),
                      pd.DataFrame(V2_RMSES)], keys=['old', 'new']).to_csv(statistics_dir+savename.split('.')[0]+'_RMSE.csv', decimal=',')

    if v1:
        axs[1].plot(np.mean(V1_RMSES.values, axis=1), color=color, linewidth=1.5, linestyle=(0, (4, 4)))
        axs[0].plot(np.mean(V1_BIASES.values, axis=1), color=color, linewidth=1.5, linestyle=(0, (4, 4)))

    if v2:
        axs[1].plot(np.mean(V2_RMSES.values, axis=1), color=color, linewidth=1.5)
        axs[0].plot(np.mean(V2_BIASES.values, axis=1), color=color, linewidth=1.5)

    axs[1].set_xlim(-1, 12)
    axs[1].set_xticks(np.linspace(0, 10, 6), month_names[::2])
    
    axs[0].set_xlim(-1, 12)
    axs[0].set_xticks(np.linspace(0, 10, 6), month_names[::2])

    if v1:
        axs[0].scatter(-0.75, np.nanmean(V1_BIASES.values), color=color, marker='v')
        axs[1].scatter(-0.75, np.nanmean(V1_RMSES.values), color=color, marker='v')
    if v2:
        axs[1].scatter(-0.25, np.nanmean(V2_RMSES.values), color=color, marker='X')
        axs[0].scatter(-0.25, np.nanmean(V2_BIASES.values), color=color, marker='X')

axs[0].set_ylabel('Bias')
axs[1].set_ylabel('RMSE')
axs[2].set_ylabel('N [% of max]')

axs[0].hlines(0, 0, 12, color='black', zorder=-1, linewidth=0.5, linestyle=(0, (1, 1)))
axs[1].hlines(0, 0, 12, color='black', zorder=-1, linewidth=0.5, linestyle=(0, (1, 1)))

legend_signs = [Line2D([0], [0], color='black', lw=4),
                Line2D([0], [0], color='#A0B1BA', lw=4),
                Line2D([0], [0], color='#00CD6C', lw=4),
                Line2D([0], [0], color='#A6761D', lw=4),
                Line2D([0], [0], marker='v', markersize=10, markerfacecolor='black', color='w'),
                Line2D([0], [0], marker='X', markersize=10, markerfacecolor='black', color='w'),
                Line2D([0], [0], color='black', linestyle=(0, (1, 1)), lw=4),
                Line2D([0], [0], color='black', lw=4)
                ]

if v1 and v2:
    axs[0].legend(legend_signs, ['Full snow cover only',
                                 'Constant pixels',
                                 'Forest',
                                 'Tundra',
                                 'Old scheme average value',
                                 'New scheme average value',
                                 'Old albedo scheme',
                                 'New albedo scheme'], prop={'size': 9})
elif v1:
    axs[0].legend(legend_signs[0:5], ['Full snow cover only',
                                 'Constant pixels',
                                 'Forest',
                                 'Tundra',
                                 'Average value',
                                 'New scheme average value',
                                 'Old albedo scheme',
                                 'New albedo scheme'], prop={'size': 12})

letters = \
    ['a', 'b', 'c']

for x in range(3):
    axs[x].text(-0.1, 1.1, letters[x], transform=axs[x].transAxes, size=10, weight='bold')

plt.savefig(fig_save_dir+savename, dpi=300)

