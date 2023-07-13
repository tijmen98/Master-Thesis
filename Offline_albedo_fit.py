from scipy.optimize import minimize
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

fig, axs = plt.subplots(3, 1, sharex=True, figsize=(10, 8))

for i, year in enumerate(['2005', '2006', '2007']):

    WEATHER = pd.read_csv('/Volumes/Tijmen/Master-Thesis/Data/Sodankyla/SODANKYLA_WEATHER_daily.csv', index_col=0).loc[
              year + '-01-01':year + '-12-31']
    RADIATION = pd.read_csv('/Volumes/Tijmen/Master-Thesis/Data/Sodankyla/' + year + '/SODANKYLA_RADIATION_daily.csv',
                            index_col=0)

    RADIATION[abs(RADIATION) < 5] = 0

    AWS_ALBEDO = abs(RADIATION['Reflected radiation (W/m2)']) / abs(
        RADIATION['Global radiation (W/m2)'])

    AWS_ALBEDO[AWS_ALBEDO > 1] = np.nan
    AWS_ALBEDO[AWS_ALBEDO < 0] = np.nan

    AWS_ALBEDO = AWS_ALBEDO.values


    def albedo_Scheme(snowfall, temperature, snowdepth, old, start_alb=float(0.75), cold_relax=0.008, warm_relax=0.24, min_alb=0.5,
                      max_alb=0.85):
        timerange = range(len(snowfall.index))
        albedo = np.zeros(len(timerange))
        albedo[0] = start_alb
        for time in timerange[0:-1]:
            if temperature[time+1] < 0:
                if old:
                    albedo[time+1] = albedo[time] - cold_relax
                if not old:
                    albedo[time + 1] = (albedo[time] - min_alb) * np.exp(-cold_relax) + min_alb
            else:
                albedo[time + 1] = (albedo[time] - min_alb) * np.exp(-warm_relax) + min_alb

            if temperature[time+1] < 2:
                albedo[time+1] = albedo[time+1] + (np.min([np.max(snowfall[time+1], 0)/10, 1]) * (max_alb - albedo[time+1]))

            if albedo[time+1] < min_alb:
                albedo[time+1] = min_alb

            if snowdepth[time+1] < 4:
                albedo[time+1] = 0.15

        return(albedo[0:-1])


    def objective(cold_relax, measured_albedo, snowfall, temperature, snowdepth, old, start_alb, warm_relax, min_alb,
                  max_alb):
        # Update the albedo_Scheme function with the tuned `cold_relax` value
        albedo = albedo_Scheme(snowfall, temperature, snowdepth, old, start_alb, cold_relax, warm_relax, min_alb, max_alb)

        # Calculate the mean squared error between measured albedo and calculated albedo
        mse = np.mean((albedo - measured_albedo) ** 2)
        return mse

    # Set the initial guess for cold_relax
    initial_guess = 0.12

    # Define the measured albedo values
    measured_albedo = AWS_ALBEDO

    # Define the other parameters required by albedo_Scheme
    snowfall = WEATHER['Precipitation amount (mm)'].astype(float)  # Define the snowfall values
    temperature = WEATHER['Air temperature (degC)'].astype(float)  # Define the temperature values
    snowdepth = WEATHER['Snow depth (cm)']  # Define the snowdepth values
    old = False  # Define the value of `old`
    start_alb = float(0.75)  # Define the value of `start_alb`
    warm_relax = 0.24  # Define the value of `warm_relax`
    min_alb = 0.5  # Define the value of `min_alb`
    max_alb = 0.85  # Define the value of `max_alb`
    cold_relax = 0.08
    figrange=[50, 100]
    fitrange=[80, 120]

    # Perform the optimization
    result = minimize(objective, initial_guess, args=(measured_albedo[fitrange[0]: fitrange[1]-1], snowfall[fitrange[0]: fitrange[1]], temperature[fitrange[0]: fitrange[1]], snowdepth[fitrange[0]: fitrange[1]], old, start_alb, warm_relax, min_alb, max_alb))

    # Retrieve the optimized value of cold_relax
    optimal_cold_relax = result.x[0]

    new_albedo = albedo_Scheme(snowfall, temperature, snowdepth, False, start_alb, 0.12, warm_relax, min_alb, max_alb)
    old_albedo = albedo_Scheme(snowfall, temperature, snowdepth, True, start_alb, 0.008, warm_relax, min_alb, max_alb)

    newRMSE = np.sqrt(np.nanmean((measured_albedo[figrange[0]: figrange[1]]-new_albedo[figrange[0]: figrange[1]])**2))
    oldRMSE = np.sqrt(np.nanmean((measured_albedo[figrange[0]: figrange[1]]-old_albedo[figrange[0]: figrange[1]])**2))

    newBIAS = np.nanmean(new_albedo[figrange[0]: figrange[1]])-np.nanmean(measured_albedo[figrange[0]: figrange[1]])
    oldBIAS = np.nanmean(old_albedo[figrange[0]: figrange[1]])-np.nanmean(measured_albedo[figrange[0]: figrange[1]])

    axs[i].plot(albedo_Scheme(snowfall, temperature, snowdepth, False, start_alb, 0.12, warm_relax, min_alb, max_alb), label='New scheme', color='red')
    axs[i].plot(albedo_Scheme(snowfall, temperature, snowdepth, True, start_alb, 0.008, warm_relax, min_alb, max_alb), label='Old scheme', color='green')
    axs[i].plot(measured_albedo, label='Measured', color='black')
    axs[i].set_xlim(figrange[0], figrange[1])
    axs[i].set_ylabel('Albedo')
    axs[i].legend()
    axs[i].set_title(year)


    print(year+' & ' + str(np.round(oldBIAS, 3))+' & '+str(np.round(newBIAS, 3))+' & '+str(np.round(abs(oldBIAS)-abs(newBIAS), 3))+' & '+str(np.round(oldRMSE, 3))+' & '+str(np.round(newRMSE, 3))+' & '+str(np.round(abs(oldRMSE)-abs(newRMSE), 3)))


axs[2].set_xlabel('Day')

plt.savefig('/Users/tijmen/Desktop/Albedo_fit/albedo_sodankyla_late.png', dpi=300)



