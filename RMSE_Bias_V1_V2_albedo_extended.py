import numpy as np
import matplotlib.pyplot as plt

fig_save_dir = '/Users/tijmen/Desktop/Figures_Thesis/'

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

fig_save_dir = '/Users/tijmen/Desktop/Figures_Thesis/'
statistics_dir = '/Volumes/Tijmen/Master-Thesis/Data/Statistics/'

arctic = True

constant = False


fig, axs = plt.subplots(2, 1, dpi=130, figsize=(16, 12))


years = ['2002', '2003', '2004', '2005']

month_names = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October',
               'November', 'December']

V1_BIASES = pd.DataFrame(columns=years, index = month_names)
V2_BIASES = pd.DataFrame(columns=years, index = month_names)

V1_RMSES = pd.DataFrame(columns=years, index = month_names)
V2_RMSES = pd.DataFrame(columns=years, index = month_names)

for i, year in enumerate(years):

    year_dir = statistics_dir+year+'/Albedo/'

    if arctic:

        savename = 'BIAS_RMSE_albedo_arctic.png'

        V1_BIAS= pd.read_csv(year_dir+'BIAS__Albedo_v1.csv', index_col=0)
        V2_BIAS= pd.read_csv(year_dir+'BIAS__Albedo_v2.csv', index_col=0)
        V1_RMSE= pd.read_csv(year_dir+'RMSE__Albedo_v1.csv', index_col=0)
        V2_RMSE= pd.read_csv(year_dir+'RMSE__Albedo_v2.csv', index_col=0)

    if constant:
        savename = 'BIAS_RMSE_albedo_constant.png'

        V1_BIAS= pd.read_csv(year_dir+'BIAS__Albedo_v1_constant_area_.csv', index_col=0)
        V2_BIAS= pd.read_csv(year_dir+'BIAS__Albedo_v2_constant_area_.csv', index_col=0)
        V1_RMSE= pd.read_csv(year_dir+'RMSE__Albedo_v1_constant_area_.csv', index_col=0)
        V2_RMSE= pd.read_csv(year_dir+'RMSE__Albedo_v2_constant_area_.csv', index_col=0)


    V1_BIASES.loc[:,year] = list(V1_BIAS.iloc[:, 0])
    V2_BIASES.loc[:,year] = list(V2_BIAS.iloc[:, 0])
    V1_RMSES.loc[:,year] = list(V1_RMSE.iloc[:, 0])
    V2_RMSES.loc[:,year] = list(V2_RMSE.iloc[:, 0])

    if i == 0:
        axs[0].plot(V1_BIAS, color='green', linewidth=0.5, linestyle=(0, (1, 1)), label='Old, one year')
        axs[0].plot(V2_BIAS, color='red', linewidth=0.5, linestyle=(0, (1, 1)), label='New, one year')
    else:
        axs[0].plot(V1_BIAS, color='green', linewidth=0.5, linestyle=(0, (1, 1)))
        axs[0].plot(V2_BIAS, color='red', linewidth=0.5, linestyle=(0, (1, 1)))

    axs[1].plot(V1_RMSE, color='green', linewidth=0.5, linestyle=(0, (1, 1)))
    axs[1].plot(V2_RMSE, color='red', linewidth=0.5, linestyle=(0, (1, 1)))


pd.concat([pd.DataFrame(V1_BIASES),
                  pd.DataFrame(V2_BIASES)], keys=['old', 'new']).to_csv(statistics_dir+savename.split('.')[0]+'_BIAS.csv', decimal=',')

pd.concat([pd.DataFrame(V1_RMSES),
                  pd.DataFrame(V2_RMSES)], keys=['old', 'new']).to_csv(statistics_dir+savename.split('.')[0]+'_RMSE.csv', decimal=',')

axs[1].plot(np.mean(V1_RMSES.values, axis=1), color='green', linewidth=1.5)
axs[1].plot(np.mean(V2_RMSES.values, axis=1), color='red', linewidth=1.5)

axs[0].plot(np.mean(V1_BIASES.values, axis=1), color='green', linewidth=1.5, label='Old, multi-year mean')
axs[0].plot(np.mean(V2_BIASES.values, axis=1), color='red', linewidth=1.5, label='New,  multi-year mean')

axs[1].set_xlim(-1, 12)
axs[1].set_xticks(np.linspace(0, 10, 6), month_names[::2])

axs[0].set_xlim(-1, 12)
axs[0].set_xticks(np.linspace(0, 10, 6), month_names[::2])

axs[0].scatter(-0.5, np.nanmean(V1_BIASES.values), color='green', marker='+', label='Old, mean value')
axs[0].scatter(-0.5, np.nanmean(V2_BIASES.values), color='red', marker='+', label='New, mean value')

axs[1].scatter(-0.5, np.nanmean(V1_RMSES.values), color='green', marker='+')
axs[1].scatter(-0.5, np.nanmean(V2_RMSES.values), color='red', marker='+')

axs[0].set_ylabel('Bias [cm]')
axs[1].set_ylabel('RMSE [cm]')

axs[0].hlines(0, 0, 12, color='black', zorder=-1, linewidth=0.5, linestyle=(0, (1, 1)))
axs[1].hlines(0, 0, 12, color='black', zorder=-1, linewidth=0.5, linestyle=(0, (1, 1)))

axs[0].legend()
axs[1].legend()

plt.savefig(fig_save_dir+savename, dpi=300)

