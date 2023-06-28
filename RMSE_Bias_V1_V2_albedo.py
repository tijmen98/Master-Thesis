import numpy as np
import matplotlib.pyplot as plt

fig_save_dir = '/Users/tijmen/Desktop/Figures_Thesis/'

V1_BIAS = [
    [0.109, 0.068, 0.035, 0.002, 0.045, 0.111, 0.303, 0.376, 0.230, 0.171, 0.132, 0.106],
    [0.068, 0.044, 0.001, 0.002, 0.051, 0.114, 0.241, 0.395, 0.268, 0.175, 0.138, 0.103],
    [0.063, 0.040, 0.001, 0.003, 0.028, 0.067, 0.306, 0.436, 0.252, 0.166, 0.122, 0.086]
            ]

V2_BIAS = [
    [0.069, 0.023, -0.013, -0.039, -0.008, 0.068, 0.307, 0.349, 0.166, 0.104, 0.089, 0.061],
    [0.044, 0.007, -0.029, -0.033, 0.010, 0.061, 0.237, 0.356, 0.219, 0.125, 0.108, 0.082],
    [0.037, 0.006, -0.032, -0.034, -0.012, 0.028, 0.291, 0.432, 0.186, 0.108, 0.086, 0.064]
            ]

V1_RMSE = [
    [0.143, 0.155, 0.156, 0.134, 0.142, 0.211, 0.317, 0.399, 0.28, 0.219, 0.174, 0.138],
    [0.133, 0.147, 0.144, 0.147, 0.145, 0.199, 0.292, 0.438, 0.333, 0.236, 0.195, 0.146],
    [0.127, 0.141, 0.141, 0.138, 0.139, 0.174, 0.360, 0.456, 0.306, 0.226, 0.180, 0.138]
            ]

V2_RMSE = [
    [0.124, 0.134, 0.134, 0.141, 0.143, 0.191, 0.311, 0.355, 0.237, 0.173, 0.144, 0.188],
    [0.114, 0.127, 0.137, 0.146, 0.138, 0.167, 0.270, 0.409, 0.287, 0.192, 0.162, 0.125],
    [0.107, 0.123, 0.137, 0.141, 0.142, 0.160, 0.327, 0.444, 0.253, 0.180, 0.149, 0.119]
            ]

fig, axs = plt.subplots(2, 1, dpi=300, figsize=(8, 6))

years = ['2002', '2003', '2004']

month_names = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October',
               'November', 'December']

for i, year in enumerate(years):

    axs[0].plot(V1_BIAS[i], color='green', linewidth=0.5, linestyle=(0, (1, 1)))
    axs[0].plot(V2_BIAS[i], color='red', linewidth=0.5, linestyle=(0, (1, 1)))

    axs[1].plot(V1_RMSE[i], color='green', linewidth=0.5, linestyle=(0, (1, 1)))
    axs[1].plot(V2_RMSE[i], color='red', linewidth=0.5, linestyle=(0, (1, 1)))

axs[1].plot(np.mean(V1_RMSE, axis=0), color='green', linewidth=1.5)
axs[1].plot(np.mean(V2_RMSE, axis=0), color='red', linewidth=1.5)

axs[0].plot(np.mean(V1_BIAS, axis=0), color='green', linewidth=1.5, label='Old, monthly mean')
axs[0].plot(np.mean(V2_BIAS, axis=0), color='red', linewidth=1.5, label='New, monthly mean')

axs[1].set_xlim(-1, 12)
axs[1].set_xticks(np.linspace(0, 10, 6), month_names[::2])

axs[0].set_xlim(-1, 12)
axs[0].set_xticks(np.linspace(0, 10, 6), month_names[::2])

axs[0].scatter(-0.5, np.mean(V1_BIAS), color='green', marker='+')
axs[0].scatter(-0.5, np.mean(V2_BIAS), color='red', marker='+')

axs[1].scatter(-0.5, np.mean(V1_RMSE), color='green', marker='+', label='Old, yearly mean')
axs[1].scatter(-0.5, np.mean(V2_RMSE), color='red', marker='+', label='New, yearly mean')

axs[0].set_ylabel('Bias [albedo]')
axs[1].set_ylabel('RMSE [albedo]')

axs[0].legend()
axs[1].legend()

plt.savefig(fig_save_dir+'BIAS_RMSE_albedo.png', dpi=300)

