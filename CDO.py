cdo = Cdo()
import xarray as xr


directory = '/Users/tijmen/Documents/Tijmen/Climate_Physics/Thesis_local/Data/Snow_cover/'

ds_measure = xr.open_dataset(directory+'Measure_2001_merged.nc')
ds_racmo = xr.open_dataset(directory+'RACMO2.4_2001-01-01_2001-12-30.nc')


#cdo.remapcon()