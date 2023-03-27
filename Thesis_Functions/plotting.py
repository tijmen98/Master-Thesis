"""plot energy fluxes at certain location as a function of time"""
        
def Plot_Energy_Flux_Graph(rlat,rlon,netlongwave,netshortwave,sensibleheatflux,latentheatflux,groundheatflux=0,refreezing=0,linestyle=1,t_start=0,t_stop=1):

    netshortwave.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).plot(label='netShortWave',color='red',linestyle=linestyle)
    netlongwave.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).plot(label='netLongWave',color='blue',linestyle=linestyle)
    sensibleheatflux.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).plot(label='SensibleHeatFlux',color='green',linestyle=linestyle)
    latentheatflux.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).plot(label='LatentHeatFlux',color='yellow',linestyle=linestyle)
    
    try:
        groundheatflux.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).plot(label='GroundHeatFlux',color='orange',linestyle=linestyle)
    except: 
        print('no groundheatflux data available')
    try: 
        refreezing.isel(rlon=rlon,rlat=rlat).sel(time=slice(t_start,t_stop)).plot(label='Refreezing',color='purple',linestyle=linestyle)
    except:
        print('no refreezing data available')
    
    

"""plot energy fluxes at a certain moment as map"""
    
def Plot_Energy_Flux_Map(time,netlongwave,netshortwave,sensibleheatflux,latentheatflux,groundheatflux,cmap='bwr',mask=1):
    
    import matplotlib.colors as colors
    import matplotlib.pyplot as plt
    
    fig , axs = plt.subplots(2,2,figsize=(25, 15))
    


    (mask*netshortwave).isel(time=time).plot(ax=axs[0,0],cmap=cmap,norm=colors.CenteredNorm())
    (mask*netlongwave).isel(time=time).plot(ax=axs[0,1],cmap=cmap,norm=colors.CenteredNorm())
    (mask*sensibleheatflux).isel(time=time).plot(ax=axs[1,0],cmap=cmap,norm=colors.CenteredNorm())
    (mask*latentheatflux).isel(time=time).plot(ax=axs[1,1],cmap=cmap,norm=colors.CenteredNorm())
    
    axs[0,0].set_title('Shortwave')
    axs[0,1].set_title('Longwave')
    axs[1,0].set_title('Sensible Heatflux')
    axs[1,1].set_title('Latent Heatflux')
    
"""plot energy fluxes at a certain moment as map with percentage colorbar limits"""
    
def Plot_Energy_Flux_Map_percent_cbar(time,netlongwave,netshortwave,sensibleheatflux,latentheatflux,groundheatflux,cmap='bwr',mask=1):
    
    import matplotlib.pyplot as plt
    
    fig , axs = plt.subplots(2,2,figsize=(25, 15))
    


    (mask*netshortwave).isel(time=time).plot(ax=axs[0,0],cmap=cmap,vmin=-2,vmax=2)
    (mask*netlongwave).isel(time=time).plot(ax=axs[0,1],cmap=cmap,vmin=-2,vmax=2)
    (mask*sensibleheatflux).isel(time=time).plot(ax=axs[1,0],cmap=cmap,vmin=-10,vmax=10)
    (mask*latentheatflux).isel(time=time).plot(ax=axs[1,1],cmap=cmap,vmin=-10,vmax=10)
    
"""plot one variable at a certain moment as map and compare"""
    
def Plot_Comparison_Map_Percentage(RACMO3,RACMO4,cmap='bwr',mask=1):
    
    import matplotlib.colors as colors
    import matplotlib.pyplot as plt
    
    fig , axs = plt.subplots(1,3,figsize=(25, 5))
    
    

    (mask*RACMO3).plot(ax=axs[0],cmap=cmap,norm=colors.CenteredNorm())
    (mask*RACMO4).plot(ax=axs[1],cmap=cmap,norm=colors.CenteredNorm())
    (mask*((RACMO4-RACMO3)/RACMO3)).plot(ax=axs[2],cmap=cmap,vmin=-2,vmax=2)
    
    axs[0].set_title('RACMO 3')
    axs[1].set_title('RACMO 4')
    axs[2].set_title('Difference [%][RC4-RC3/RC3]')   
    title_time=str(RACMO4.isel(time=4).time.values)
    fig.suptitle(RACMO4.long_name+' at '+title_time.split(sep='T')[0])
"""plot one variable at a certain moment as map and compare"""
    
def Plot_Comparison_Map_Absolute(RACMO3,RACMO4,cmap='bwr',mask=1):
    
    import matplotlib.colors as colors
    import matplotlib.pyplot as plt
    
    fig , axs = plt.subplots(1,3,figsize=(25, 5))
    
    

    (mask*RACMO3).plot(ax=axs[0],cmap=cmap,norm=colors.CenteredNorm())
    (mask*RACMO4).plot(ax=axs[1],cmap=cmap,norm=colors.CenteredNorm())
    (mask*RACMO4-RACMO3).plot(ax=axs[2],cmap=cmap,vmin=-2,vmax=2)

    axs[0].set_title('RACMO 3')
    axs[1].set_title('RACMO 4')
    axs[2].set_title('Absolute difference [RC4-RC3]')   
    title_time=str(RACMO4.isel(time=4).time.values)
    fig.suptitle(RACMO4.long_name+' at '+title_time.split(sep='T')[0])