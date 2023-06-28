import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

fig_save_dir = '/Users/tijmen/Desktop/Figures_Thesis/'
statistics_dir = '/Volumes/Tijmen/Master-Thesis/Data/Statistics/'

arctic = False

norway = False

siberia = True

V1_BIAS = [
    [25.60, 31.19, 33.63, 28.15, 12.59, 3.85, 1.39, 0.91, 1.14, 4.56, 12.00, 17.54],
    [25.446, 29.655, 27.818, 24.833, 14.331, 4.165, 1.381, 0.900, 1.365, 3.374, 8.732, 17.550],
    [23.827, 28.113, 29.186, 26.324, 13.268, 4.295, 1.310, 0.812, 1.327, 3.405, 10.439, 19.505]
            ]

V2_BIAS = [
    [23.31, 30.93, 32.90, 26.72, 13.35, 3.45, 1.31, 0.85, 1.09, 4.54, 11.89, 17.44],
    [25.235, 29.506, 27.399, 23.629, 12.970, 3.639, 1.274, 0.837, 1.302, 3.293, 8.661, 17.535],
    [23.623, 25.581, 29.213, 25.424, 12.431, 3.928, 1.226, 0.784, 1.301, 3.303, 10.330, 19.396]
            ]

V1_RMSE = [
    [51, 63.7, 68.6, 63.7, 49.1, 28.4, 12.1, 10, 7.6, 13.8, 22.7, 32.4],
    [47.500, 55.700, 56.200, 55.800, 44.000, 24.000, 11.700, 9.200, 6.900, 14.700, 19.000, 35.800],
    [46.3, 54.300, 59.300, 61.200, 44.200, 24.200, 11.000, 6.700, 6.300, 10.000, 24.100, 39.600]
            ]

V2_RMSE = [
    [50.50, 63.20, 67.90, 61.90, 46.90, 27.10, 11.10, 9.60, 7.30, 13.70, 22.40, 32.00],
    [47.00, 55.20, 55.60, 54.00, 41.40, 21.90, 10.50, 8.60, 6.40, 14.40, 18.70, 35.60],
    [46.00, 54.70, 59.30, 60.10, 42.90, 22.90, 10.20, 6.70, 6.00, 9.80, 24.00, 39.50]
            ]

fig, axs = plt.subplots(2, 1, dpi=300, figsize=(8, 6))

years = ['2002', '2003', '2004']

month_names = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October',
               'November', 'December']

V1_BIASES = pd.DataFrame(columns=years, index = month_names)
V2_BIASES = pd.DataFrame(columns=years, index = month_names)

V1_RMSES = pd.DataFrame(columns=years, index = month_names)
V2_RMSES = pd.DataFrame(columns=years, index = month_names)

for i, year in enumerate(years):

    year_dir = statistics_dir+year+'/Snowdepth/'

    if arctic:

        savename = 'BIAS_RMSE_snowdepth_arctic.png'

        V1_BIAS= pd.read_csv(year_dir+'BIAS_arctic_snowdepth_v1.csv', index_col=0)
        V2_BIAS= pd.read_csv(year_dir+'BIAS_arctic_snowdepth_v2.csv', index_col=0)
        V1_RMSE= pd.read_csv(year_dir+'RMSE_arctic_snowdepth_v1.csv', index_col=0)
        V2_RMSE= pd.read_csv(year_dir+'RMSE_arctic_snowdepth_v2.csv', index_col=0)

    if norway:
        savename = 'BIAS_RMSE_snowdepth_norway.png'

        V1_BIAS= pd.read_csv(year_dir+'BIAS_norway_snowdepth_v1.csv', index_col=0)
        V2_BIAS= pd.read_csv(year_dir+'BIAS_norway_snowdepth_v2.csv', index_col=0)
        V1_RMSE= pd.read_csv(year_dir+'RMSE_norway_snowdepth_v1.csv', index_col=0)
        V2_RMSE= pd.read_csv(year_dir+'RMSE_norway_snowdepth_v2.csv', index_col=0)

    if siberia:
        savename = 'BIAS_RMSE_snowdepth_siberia.png'

        V1_BIAS= pd.read_csv(year_dir+'BIAS_syberia_snowdepth_v1.csv', index_col=0)
        V2_BIAS= pd.read_csv(year_dir+'BIAS_syberia_snowdepth_v2.csv', index_col=0)
        V1_RMSE= pd.read_csv(year_dir+'RMSE_syberia_snowdepth_v1.csv', index_col=0)
        V2_RMSE= pd.read_csv(year_dir+'RMSE_syberia_snowdepth_v2.csv', index_col=0)


    V1_BIASES.loc[:,year] = list(V1_BIAS.iloc[:, 0])
    V2_BIASES.loc[:,year] = list(V2_BIAS.iloc[:, 0])
    V1_RMSES.loc[:,year] = list(V1_RMSE.iloc[:, 0])
    V2_RMSES.loc[:,year] = list(V2_RMSE.iloc[:, 0])

    axs[0].plot(V1_BIAS, color='green', linewidth=0.5, linestyle=(0, (1, 1)))
    axs[0].plot(V2_BIAS, color='orange', linewidth=0.5, linestyle=(0, (1, 1)))

    axs[1].plot(V1_RMSE, color='green', linewidth=0.5, linestyle=(0, (1, 1)))
    axs[1].plot(V2_RMSE, color='orange', linewidth=0.5, linestyle=(0, (1, 1)))

axs[1].plot(np.mean(V1_RMSES.values, axis=1), color='green', linewidth=1.5)
axs[1].plot(np.mean(V2_RMSES.values, axis=1), color='orange', linewidth=1.5)

axs[0].plot(np.mean(V1_BIASES.values, axis=1), color='green', linewidth=1.5, label='Version 1 monthly mean')
axs[0].plot(np.mean(V2_BIASES.values, axis=1), color='orange', linewidth=1.5, label='Version 2 monthly mean')

axs[1].set_xlim(-1, 12)
axs[1].set_xticks(np.linspace(0, 10, 6), month_names[::2])

axs[0].set_xlim(-1, 12)
axs[0].set_xticks(np.linspace(0, 10, 6), month_names[::2])

axs[0].scatter(-0.5, np.mean(V1_BIASES.values), color='green', marker='+')
axs[0].scatter(-0.5, np.mean(V2_BIASES.values), color='orange', marker='+')

axs[1].scatter(-0.5, np.mean(V1_RMSES.values), color='green', marker='+', label='Version 1 yearly mean')
axs[1].scatter(-0.5, np.mean(V2_RMSES.values), color='orange', marker='+', label='Version 2 yearly mean')

axs[0].set_ylabel('Bias')
axs[1].set_ylabel('RMSE')

if siberia:
    axs[0].set_ylim(-10, 10)
    axs[1].set_ylim(0, 30)

axs[0].legend()
axs[1].legend()

plt.savefig(fig_save_dir+savename, dpi=300)
