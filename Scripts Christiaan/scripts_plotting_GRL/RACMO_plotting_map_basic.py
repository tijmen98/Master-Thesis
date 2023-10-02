# -*- coding: utf-8 -*-
"""
Script to plot a map of Greenland for any variable for RACMO and a reference RACMO run.

CREATED BY: Christiaan van Dalum

"""
### Loading packages
import numpy as np
import os
from mpl_toolkits.basemap import Basemap
import xarray as xr
import pandas as pd
import sys
sys.path.insert(1, '/Users/dalum/Documents/Scripts/functions/')
from grlmap_basic import grlmap
from get_latlon_shifted import get_latlon_shifted
from own_colormap import own_colormap


### Names and directories
savedir_base            = "/Users/dalum/Resultaten/RACMO24/FGRN11/RACMO24_1_complete2_md_debug/" #Save directory
RACMOdir_base           = "/Users/dalum/RACMO2.4/experiments/RACMO24_1_complete2_md_debug/FGRN11/NC_MD/"  # Name of the directory with the RACMO data
RACMOfilename_base      = "FGRN11.RACMO24_1_complete2_md_debug.DD" # Base name of the files
icemask_loc             = "/Users/dalum/clim_files/RACMO2.4_climfields/FGRN11/netcdf/" # Location of the icemask, use the climate fields for this
icemask_name            = "cl00090000"
Rglacmin                = 0.5       # Minimum ice tile fraction. Set to 0.5 for now

### With a reference run, e.g. RACMO2.3p3
ref                    = 1     # Switch to turn on the reference run.  1 for ON                                                                                     # Reference run on, e.g. RACMO2.3p2
RACMOrefdir_base       = "/Users/dalum/RACMO_output/RACMO2.3p3/FGRN11/NC_MD/"            # Name of the directory of the reference  (ONLY WITH REF = 1)
RACMOreffilename_base  = 'FGRN11.CVD_SVN548_bia_update_FINAL_4_imp5_NEW.DD'          # Base name of the files of the reference  (ONLY WITH REF = 1)

### Time and gridpoint settings
start_time              = '2000-09-01'                   # Start time 'YYYY-MM-DD'
end_time                = '2000-12-31'                   # End time 'YYYY-MM-DD'
time_forloop            = [2000] #[2001,2011] #[2001,2011]   # Time elements that we have to go through to open the correct data files         

### Other settings
var_name                = ['smbgl']  # Variable name      ##
# For RACMO2.4:
#['smbgl','pr','sublgl','sndiv','totrunoff','mltgl','refreeze']    #['smbgl','pr','sublgl','sndiv','totrunoff','mltgl','refreeze']                       
var_name_ref            = ['smb']    # Variable name      ##

mm_we_year              = 0     # Switch if we want to convert the units to mm we yr-1
mm_we_months            = 1     # 1 for mm w.e. per month.  
nummonths               = 3     # Number of months.


### The bounds for the colors:
# for the normal figures:
#bounds = np.array([-2000,-1750,-1500,-1250,-1000,-750,-500,-400,-300,-200,-100,0,100,200,300,400,500,750,1000,1250,1500,1750,2000])  
#bounds                    = np.array([0,1,5,10,30,50,100,200,300,400])    
bounds    = np.array([-500, -400, -300, -200, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 200, 300, 400, 500])    # Steps for the colormap of the figure
cbarticks = [-500,-100,0,100,500]  #[-2000,-500,0,500,2000]  #[-4000,-1000,0,1000,4000] #[-2000,-500,0,500,2000]     # The ticks to plot
#cbarticks = bounds

# for the diff figures
#bounds_diff = np.array([-500, -450, -400, -350, -300, -250, -200, -150, -125, -100, -75, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500])
bounds_diff = np.array([-500, -400, -300, -200, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 200, 300, 400, 500])   # Steps for the colormap of the figure
cbarticks_diff = [-500,-100,0,100,500]   # The ticks to plot

# for the percentage figures
bounds_percent = np.array([-100,-90,-80,-70,-60, -50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100])   # Steps for the colormap of the figure
cbarticks_percent = [-100,-50,0,50,100]   # The ticks to plot

# For RACMO2.3p3:
#['alb', 'ci', 'clcov', 'clcovc', 'clcovH', 'clcovL', 'clcovM', 'evap', 'evapl', 'ff10m', 'ff10max', 'ffgmax','frsnr', 'gbot',
# 'hgtsrf', 'latf', 'lw0n', 'lw0ncf', 'lwsd', 'lwsn', 'lwsncf', 'lz0h', 'meltin', 'mslp', 'mu0', 'preccv', 'precip', 'precls',
# 'psurf', 'q2m', 'qci', 'qii', 'qli', 'qvi', 'raind', 'refreeze', 'rh2m', 'runoff', 'senf', 'smb', 'sndiv', 'snowcv', 'snowfall',
# 'snowmass', 'snowmelt', 'snowrho', 'snowtemp', 'sst', 'subl', 'suds', 'sund', 'sus', 'sw0d', 'sw0n', 'sw0ncf', 'sw0u', 'swsd',
# 'swsn', 'swsncf', 'swsu', 't2m', 't2max', 't2min', 'taucvi', 'tileFR3', 'tileFR4', 'tileFR5', 'tileFR6', 'tileFR7', 'tileFR8',
# 'totpore', 'totwat', 'trds', 'trds_rsn', 'trds_x', 'trds_y', 'tsal1', 'tsal2', 'tskin', 'tsmax', 'tsmin', 'turbdis', 'u10m',
# 'ustr', 'ustrgw', 'v10m', 'vdisgw', 'vstr', 'vstrgw', 'z0m', 'zsnow']
#
## SMB components: smb, precip, subl, sndiv, runoff, snowmelt, refreeze.


timeyear = 3600*24*365.25       # Time in seconds of a year
waterdensity = 1000.            # Water density

#%%

#### Making the Save directory
savedir = savedir_base + 'maps/'
if not os.path.exists(savedir):                                                     # Making the save directory if it does not exist
    os.makedirs(savedir)   

### Define time as a pandas time series
time_series = pd.date_range(start=start_time, end=end_time, freq='1D')
time_series = xr.Dataset({'time': ('time', time_series)})


# Load the ice mask and do some corrections
icemask = xr.open_dataset("%s/%s.nc"%(icemask_loc,icemask_name))
icemask = np.array(icemask['var91'][0,0])
# Set the icemask to 1 and 0.
icemask2 = np.where(icemask > Rglacmin, 1.0, icemask)    
icemask2 = np.where(icemask2 < Rglacmin, 0.0, icemask2)  
icemask3 = np.where(np.isnan(icemask2), 0.0, icemask2) 


def extract_data(RACMOdir, RACMO_filename, varname,time_series, atts = 1):
    """
    Function to extract the data from RACMO netcdf files
    
    ### INPUT ###
    RACMOdir: Directory of the RACMO data.
    RACMO_filename: RACMO file name
    varname: Variable name
    time_series: pandas time series for the time to consider
    atts: switch to extract attributes: the savename, ylabel and unit (DEFAULT: 1)
    
    ### OUTPUT ###
    data: xarray dataset with the considered variable
    savename: Standard name to use for saving the figure
    ylabel: Name for the ylabel
    units: Units for the ylabel
    
    """

    ### Open the data set and extract the variable and grid point
    data = xr.open_dataset(RACMOdir + varname + '%s.nc'%(RACMO_filename))[varname][:,0,:,:]
    ### Select the correct time frame
    data = data.sel(time=slice(time_series["time"][0],time_series["time"][-1]))    
    ### Extract the attributes
    if atts == 1:
        savename = data.attrs['standard_name']
        ylabel = data.attrs['long_name']
        units = data.attrs['units']
    ### Extent the data set and fill missing values
    if atts == 1:
        return data,savename,ylabel,units
    else:
        return data
 
all_data_ref = []; all_diff = []; all_diff_percent = []
all_domain_averaged_percent = []; all_data = []; all_units = []
  
for j in range(0,len(var_name)):
    data = np.atleast_1d(0)                # Initialize
    
    ### Read the data and manipulate it a bit.
    for i in range(0,len(time_forloop)):
        if os.path.isfile(RACMOdir_base + var_name[j] + '.KNMI-%s.%s.nc'%(time_forloop[i],RACMOfilename_base)) == False:
            print ("%s does not exist... exit"%(var_name[j] + '.KNMI-%s.%s.nc'%(time_forloop[i],RACMOfilename_base)))
            break
        else:
            data_temp, savename, ylabel, units = extract_data(RACMOdir_base,'.KNMI-%s.%s'%(time_forloop[i],RACMOfilename_base),var_name[j],time_series) 
            # Save lat/lon data on the first time step       
            if i == 0  and j==0:
                lat = data_temp["lat"]; lon = data_temp["lon"]  
            rlat_cor = np.round(data_temp.coords["rlat"].values,decimals=1)
            rlon_cor = np.round(data_temp.coords["rlon"].values,decimals=2)
            data_temp = data_temp.drop({"lon","lat","height"})
            data_temp = data_temp.assign_coords(rlat=rlat_cor, rlon=rlon_cor)
            if len(data) != 1:
                data = xr.merge([data,data_temp])
            else:
                data = data_temp
        ### Extent the data set and fill missing values
        data = xr.merge([data,time_series]) 
        data = data[var_name[j]]*icemask2
    
    # Convert the data to mm w.e. yr-1 if necessary
    if mm_we_year == 1:
        data = data/waterdensity * timeyear *  1000
        
    data_mean = data.mean(dim='time')
        
    # Convert the data to mm w.e. month-1 if necessary
    if mm_we_months == 1:
        data = data.sum(dim='time')
        data = data/waterdensity * 1000 *  3600*24
        data_mean = data/nummonths
        
    all_data.append(data_mean)
    all_units.append(units)
    
    ### Read the reference data and manipulate it a bit.
    
    if ref == 1:
        data_ref = np.atleast_1d(0)           # Initialize
        for i in range(0,len(time_forloop)):
            if os.path.isfile(RACMOrefdir_base + var_name_ref[j] + '.KNMI-%s.%s.nc'%(time_forloop[i],RACMOreffilename_base)) == False:
                print ("%s does not exist... exit"%(var_name_ref[j] + '.KNMI-%s.%s.nc'%(time_forloop[i],RACMOreffilename_base)))
                pass
            else:
                data_ref_temp = extract_data(RACMOrefdir_base,'.KNMI-%s.%s'%(time_forloop[i],RACMOreffilename_base),var_name_ref[j],time_series,atts = 0)  
                rlat_cor = np.round(data_ref_temp.coords["rlat"].values,decimals=1)
                rlon_cor = np.round(data_ref_temp.coords["rlon"].values,decimals=2)
                data_ref_temp = data_ref_temp.drop({"lon","lat","height"})
                data_ref_temp = data_ref_temp.assign_coords(rlat=rlat_cor, rlon=rlon_cor)
                if len(data_ref) != 1:
                    data_ref = xr.merge([data_ref,data_ref_temp])
                else:
                    data_ref = data_ref_temp
                    
        data_ref = xr.merge([data_ref,data_ref_temp])
        data_ref = data_ref[var_name_ref[j]] * icemask2
        
        
        if mm_we_year == 1:
            data_ref = data_ref/waterdensity * timeyear *  1000
        data_ref_mean = data_ref.mean(dim='time')
        
        
        if mm_we_months == 1:
            data_ref = data_ref.sum(dim='time')
            data_ref = data_ref/waterdensity * 1000 *  3600*24
            data_ref_mean = data_ref/nummonths
            
        diff = data_mean - data_ref_mean
        diff_percent = diff / abs(data_ref_mean) * 100
        domain_averaged_percent = np.nansum(diff.values)/abs(np.nansum(data_ref_mean.values))*100
    
        all_data_ref.append(data_ref_mean)
        all_diff.append(diff)
        all_diff_percent.append(diff_percent)
        all_domain_averaged_percent.append(domain_averaged_percent)
            
### Do some minor corrections to the lat/lon data.
block2 = np.zeros(40)
block2[:] = (50, 0, 0, 0, 0, 10, 326, 0, 327, 0, -15800, 0, 0, -17650, 0, 0, 136, 
16800, 0, 0, 14850, 0, 0, 180, 0, 180, 0, 64, 0, 0, 0, 0, -18000, 0, 0, 
-37500, 0, 0, 1, 1 )
gridlatlon = get_latlon_shifted(block2)
lat2 = np.zeros_like(lat)    
lon2 = np.zeros_like(lon)
for i in range(0,len(lat[0,:])):
    lat2[:,i] = gridlatlon[0,i,:]
    lon2[:,i] = gridlatlon[1,i,:]    
   
### Making the base map
m = Basemap(lat_0 = 57,lon_0=-42, llcrnrlon=-55, llcrnrlat=57, urcrnrlon=14, urcrnrlat=82, resolution='h',projection='stere')  
x,y = m(lon2,lat2)  


#%% 
###
### PLOTTING
####
### Map of the variable    

for j in range(0,len(var_name)):
    var_name_plot = var_name[j]
    # var_name_plot           = ['smb','precipitation','sublimation','sndiv','totrunoff','mltgl','refreeze'] 
    savename = var_name[j]
    ### Settings
    if mm_we_year == 1:
        units = 'mm w.e. yr$^{-1}$'      # Units
        steps = len(bounds)-1            # Number of steps for the color map 
        lowval = -2000                   # Lowest value to plot
        highval = 2000                   # Highest value to plot
        cbarticks_local = cbarticks      # the ticks of the colorbars
        cmap = own_colormap(totalcolors=len(bounds)-1)      # The color map
        bounds_local = bounds            # The bounds of the color map
    elif  mm_we_months == 1: 
        units = 'mm w.e. month$^{-1}$'
        steps = len(bounds)-1
        lowval = -500
        highval = 500
        cbarticks_local = cbarticks
        cmap = own_colormap(totalcolors=len(bounds)-1)   
        bounds_local = bounds
    else:
        units = all_units[j]
        steps = 18
        lowval = np.nanmin(all_data[j])
        highval = np.nanmax(all_data[j])
        cbarticks_local = 4
        cmap = own_colormap(totalcolors=20)
        bounds_local = ''

    # Plot the map
    grlmap(all_data[j],m,x,y,savedir, savename + '%s-%s'%(start_time[0:4],end_time[0:4]), lowval = lowval, highval = highval, 
           cbarticks=cbarticks_local, titlefontsize = 18,cbartitle = var_name_plot + " [%s]"%(units),cbarticksize = 15,
           colormap = cmap,steps=steps, bounds = bounds_local,meanon = True,
           mean = np.round(np.nanmean(all_data[j]),decimals=1), icemask = icemask3, contourcolor = 'black')  
        
    ### Plot the reference maps
    if ref == 1:
        ### Map of the variable -- REFERENCE
        savedir2 = savedir + 'RACMOref/'
        if not os.path.exists(savedir2):                                                     # Making the save directory if it does not exist
            os.makedirs(savedir2)   
        
        savename = var_name_ref[j]
        ### Settings
        if mm_we_year != 1 and mm_we_months != 1:
            lowval = np.nanmin(all_data_ref[j])
            highval = np.nanmax(all_data_ref[j])
            
        grlmap(all_data_ref[j],m,x,y,savedir2, savename + '%s-%s'%(start_time[0:4],end_time[0:4]), lowval = lowval, highval = highval, 
               cbarticks=cbarticks_local, titlefontsize = 18,cbartitle = var_name_plot + " [%s]"%(units),cbarticksize = 15,
               colormap = cmap,steps=steps, bounds = bounds_local,meanon = True,
               mean = np.round(np.nanmean(all_data_ref[j]),decimals=1), icemask = icemask3, contourcolor = 'black')   
        
    ### Map of the the difference: current - reference
        savedir3 = savedir + 'diff/'
        if not os.path.exists(savedir3):                                                     # Making the save directory if it does not exist
            os.makedirs(savedir3)   
        
        ### Load own colormap for difference plots
        ### Settings
        if mm_we_year == 1:
            steps = len(bounds_diff)-1
            lowval = -500
            highval = 500
            cbarticks_local = cbarticks_diff
            cmap = own_colormap(totalcolors=len(bounds_diff)-1) 
            bounds_local = bounds_diff
        elif  mm_we_months == 1: 
            steps = len(bounds_diff)-1
            lowval = -500
            highval = 500
            cbarticks_local = cbarticks_diff
            cmap = own_colormap(totalcolors=len(bounds_diff)-1) 
            bounds_local = bounds_diff
        else:
            steps = 18
            lowval = np.nanmin(all_diff[j])
            highval = np.nanmax(all_diff[j])
            cbarticks_local = 4
            cmap = own_colormap(totalcolors=20)
            bounds_local = ''
            
        # Plot the diff
        grlmap(all_diff[j],m,x,y,savedir3, 'diff' + savename + '%s-%s'%(start_time[0:4],end_time[0:4]), lowval = lowval, highval = highval,
               cbarticks=cbarticks_local, titlefontsize = 18,cbartitle = "$\Delta$%s [%s]"%(var_name_plot,units), cbarticksize = 15,
               colormap = cmap,meanon = True, mean = np.round(np.nanmean(all_diff[j]),decimals=1),
               steps=steps, bounds = bounds_local, icemask = icemask3, contourcolor = 'black')   
        
        ### Settings
        if mm_we_year == 1:
            steps = len(bounds_percent)-1
            lowval = -100
            highval = 100
            cbarticks_local = cbarticks_percent
            cmap = own_colormap(totalcolors=len(bounds_percent)-1,paper2=True)  
            bounds_local = bounds_percent
        elif  mm_we_months == 1: 
            steps = len(bounds_percent)-1
            lowval = -100
            highval = 100
            cbarticks_local = cbarticks_percent
            cmap = own_colormap(totalcolors=len(bounds_percent)-1,paper2=True)  
            bounds_local = bounds_percent
        else:
            steps = 18
            lowval = -100
            highval = 100
            cbarticks_local = 4
            cmap = own_colormap(totalcolors=20)
            bounds_local = ''
            
        ### Diff in percentage with respect to the reference
        grlmap(all_diff_percent[j],m,x,y,savedir3, 'diff_percent' + savename + '%s-%s'%(start_time[0:4],end_time[0:4]), lowval = lowval, highval = highval,
               cbarticks=cbarticks_local, titlefontsize = 18,cbartitle = "$\Delta$%s [%%]"%(var_name_plot),cbarticksize = 15,
               colormap = cmap,meanon = True, mean = np.round(all_domain_averaged_percent[j],decimals=1),
               steps=steps,bounds=bounds_local, icemask = icemask3, contourcolor = 'black')   
    
