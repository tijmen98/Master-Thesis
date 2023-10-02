def grlmap(data,m,x,y,savedir,figurename,meanon = False, mean = "",
           fontsizemeridians = 14,lowval=0,highval=1,steps=101,colormap="RdYlBu_r",
           cbarticksize=14,cbartitle=" ",cbarticks=6,title=" ",
           titlefontsize=15, extend = 'both', ticktype = float, bounds = '',
           icemask = '',contourcolor = 'palegreen'):
    """
    grlmap: A function for making a map of greenland
    MADE BY: Christiaan van Dalum
    
    ### INPUT
     - data: Array or xarray of the data (lat,lon)
     - m: Map as is produced by basemap
     - x: x-values as is produced by basemap (lat,lon)
     - y: y-values as is produced by basemap (lat,lon)
     - savedir: Directory to save the figure in
     - figurename: Name of the figure
     
     # OPTIONAL
     - meanon: True if the mean value is plotted in the corner (DEFAULT: FALSE)
     - mean: Mean value to be plotted with meanon (DEFAULT: "")
     - fontsizemeridians: Font size of the meridians (DEFAULT: 14)
     - lowval: Lowest value to plot (DEFAULT: 0)
     - highval: Highest value to plot (DEFAULT: 1)
     - steps: Number of steps in the colorbar (DEFAULT: 101)
     - colormap: Colormap to use, can be an own color map or a standard one (DEFAULT: RdYlBu_r)
     - cbarticksize: Size of the colorbar ticks (DEFAULT: 14)
     - cbartitle: Title of the colorbar (DEFAULT: "")
     - cbarticks: Number of ticks of the colorbar (DEFAULT: 6)
     - title: Title of the plot (DEFAULT: " ")
     - titlefontsize: Font size of the title (DEFAULT: 15)
     - extend: Whether to extend the colorbar beyond the provided limits (DEFAULT: 'both')
     - ticktype: What type to use for the colorbar ticks (DEFAULT: float)
     - bounds: The values to use in the colorbar (DEFAULT: '')
     - icemask: To make a contour of the ice sheet using an ice mask (lat,lon) (DEFAULT: '')
     - contourcolor: Color to use for the ice sheet contours (DEFAULT: 'palegreen')

    ### OUTPUT
     - A figure saved as a png in the save directory
     - A figure saved as a pdf in the save directory
    """

    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib.colors as mcolors    
    
    fig,ax = plt.subplots(figsize=(10,10))
    
    # Plot meridians and coastlines
    m.drawcoastlines(color='#6D5F47', zorder = 0)
    m.drawcountries(color='#6D5F47')
    m.drawmeridians(np.arange(-60.,60,15.),labels=[0,0,0,1],fontsize=fontsizemeridians)
    m.drawparallels(np.arange(0.,90,5),labels=[1,0,0,0],fontsize=fontsizemeridians)        
                
    #Making a map, using the data and plotting it. 
    if isinstance(cbarticks,int):
        ticks = np.linspace(lowval,highval,num=cbarticks,dtype=ticktype)
    else:
        ticks = cbarticks
    
    # Colormap settings
    cmap = plt.cm.get_cmap(colormap,steps)
    
    # Plot the ice contour and the data
    if bounds == '':
        if icemask != '':
            m.contour(x,y,icemask,alpha = 0.43,levels=[0],colors = contourcolor) 
        cs_new = m.pcolormesh(x,y,data,vmin = lowval,vmax = highval, cmap = cmap)
    else:
        if icemask != '':
            m.contour(x,y,icemask,alpha = 0.43,levels=[0],colors = contourcolor)  
        norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=steps)   
        cs_new = m.pcolormesh(x,y,data,vmin = lowval,vmax = highval, cmap = cmap,norm=norm) 
        
    # Colorbar settings
    cbar = m.colorbar(cs_new, pad="5%",ticks=ticks,location="bottom", extend=extend)
    cbar.ax.tick_params(labelsize=cbarticksize,)
    cbar.set_label(cbartitle,fontsize=titlefontsize)
    
    # Title
    plt.title(title,fontsize=titlefontsize)
    
    # Plot the mean if requested
    if meanon == True:
        plt.annotate('Mean: %s'%(mean),xy=(20000,20000),fontsize = titlefontsize)

    # Save the figure    
    plt.savefig("%s/%s.png" %(savedir,figurename),bbox_inches='tight', dpi=500)
    plt.savefig("%s/%s.pdf" %(savedir,figurename),bbox_inches='tight')
    plt.close()  
    
    
    