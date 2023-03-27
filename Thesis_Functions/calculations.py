"""returns the Watt convertion for a certain cumulative joule value based on the length of the month"""

def CumJoule_to_Watt_Day(data):
    
    import xarray as xr
    
    watt = xr.full_like(data,0)
    from calendar import monthrange  
    for index in range(len(data['time'])):
        time = (data['time'][index].values.astype('str').split(sep='-')[0:2])
        totalsec = 24*60*60
        watt[index] = data.isel(time = index).values/totalsec
        watt.attrs['cell_method'] = 'time: average'
        watt.attrs['units'] = 'W m-2'
    return (watt)

def CumJoule_to_Watt_Month(data):
    
    import xarray as xr
    
    watt = xr.full_like(data,0)
    from calendar import monthrange  
    for index in range(len(data['time'])):
        time = (data['time'][index].values.astype('str').split(sep='-')[0:2])
        lenmonth = (monthrange(int(time[0]),int(time[1]))[1])
        totalsec = lenmonth*24*60*60
        watt[index] = data.isel(time = index).values/totalsec
        watt.attrs['cell_method'] = 'time: average'
        watt.attrs['units'] = 'W m-2'
    return (watt)


def Consecutive_Tundra(tile5,threshold=0.2):

    import scipy     
    import numpy as np
    import xarray as xr

    tile5_filtered = tile5.where(tile5<=threshold,1).where(tile5>=threshold)
    
    low_moments = scipy.signal.find_peaks(tile5.sum(('rlat','rlon')).transpose().to_numpy()[0])[0]
    number_of_days = xr.DataArray(np.zeros((len(tile5['rlat']),len(tile5['rlon']),int(len(low_moments)-1))))
    
    for index in range(len(low_moments)-1):
        number_of_days[:,:,index] = tile5_filtered.sel(time=slice(tile5['time'][low_moments[index]].values,tile5['time'][low_moments[index+1]].values)).sum('time').squeeze()
    
    number_of_days = number_of_days.rename(dim_0=('rlat'),dim_1=('rlon'),dim_2=('time'))
    return(number_of_days)


def return_index_from_coordinates(lat,lon,data):
    
    import numpy as np
    
    index = np.where((abs(data['lat']-lat)+abs(data['lon']-lon))==np.amin((abs(data['lat']-lat)+abs(data['lon']-lon))))

    return([int(index[1]),int(index[0])])


def Net_Radiation(net_shortwave,net_longwave,Sensible_heatflux,Latent_heatflux,Refreezing,Ground_heatflux):

    
    netradiation = net_shortwave+net_longwave+Sensible_heatflux+Latent_heatflux
    
    try:
        netradiation+Refreezing
    except:
        print('No refreezing data')

    try:
        netradiation+Ground_heatflux
    except:
        print('No ground heatflux data')
    return(netradiation)