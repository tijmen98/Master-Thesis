"""
Made by WJB on 12 Nov 2017
Information about settings to use are on ECMWF, $HOME/LBC/ToRacmo2/setttings/SETTINGS_"GRID NAME"
For FSG180:
(50, 0, 0, 0, 0, 10, 86, 0, 80, 0, -7020, 0, 0, -7560, 0, 0, 136,
    7200, 0, 0, 7740, 0, 0, 180, 0, 180, 0, 64, 0, 0, 0, 0, -26500, 0, 0,
    -44000, 0, 0, 1, 1 )
For FGRN11:
(50, 0, 0, 0, 0, 10, 326, 0, 327, 0, -15800, 0, 0, -17650, 0, 0, 136,
    16800, 0, 0, 14850, 0, 0, 180, 0, 180, 0, 64, 0, 0, 0, 0, -18000, 0, 0,
    -37500, 0, 0, 1, 1 )
For ANT27y18:
(50, 0, 0, 0, 0, 10, 262, 0, 240, 0, -30000, 0, 0, -32750, 0, 0, 136,
    29750, 0, 0, 32500, 0, 0, 180, 0, 180, 0, 64, 0, 0, 0, 0, 180000, 0, 0,
    10000, 0, 0, 1, 1 )
For ANT11:
(50, 0, 0, 0, 0, 10, 646, 0, 580, 0, -29000, 0, 0, -32300, 0, 0, 136,
    28900, 0, 0, 32200, 0, 0, 180, 0, 180, 0, 64, 0, 0, 0, 0, 180000, 0, 0,
    10000, 0, 0, 1, 1 )
"""

import numpy as np

def get_latlon_shifted(block2):
    nx    = np.int(block2[6])
    ny    = np.int(block2[8])
    # the following numbers are directly converted from thousands of a degree into
    # radians
    tdg2pi = np.pi/(1000.*180.)
    south = block2[10]*tdg2pi
    west  = block2[13]*tdg2pi
    north = block2[17]*tdg2pi
    east  = block2[20]*tdg2pi
    dx    = block2[23]*tdg2pi       # dy is equal to dx
    polat =(block2[32] + 90000.)*tdg2pi # +90 is easier for the calculations
    polon = block2[35]*tdg2pi
    
    # put these to values to zero to check if you get the real grid box centers
    latshift = -dx*0.5
    lonshift = -dx*0.5
    
    rlat   = np.zeros( (nx,ny) )
    rlon   = np.zeros( (nx,ny) )
    for ix in range(nx):
        rlat[ix,:] = np.linspace(south+latshift,north+latshift,ny)
    for iy in range(ny):
        rlon[:,iy] = np.linspace(west+lonshift,east+lonshift,nx)

# rewrite lat-lon into a unit vector (assuming round Earth, as RACMO does as well)  
    ex     = np.cos(rlat)*np.cos(rlon)
    ey     = np.cos(rlat)*np.sin(rlon)
    ez     = np.sin(rlat)

# rotate this vector for the desired projection (assuming round Earth)
# note [if you would like to figure out what is done] this are two steps into one
    gx     = (np.cos(polat)*ex - np.sin(polat)*ez)*np.cos(polon) - np.sin(polon)*ey
    gy     = (np.cos(polat)*ex - np.sin(polat)*ez)*np.sin(polon) + np.cos(polon)*ey
    gz     =  np.sin(polat)*ex + np.cos(polat)*ez
    
# convert back into real world lat lon  
    gridlatlon = np.zeros( (2,nx,ny) )  
## don't take for granted that gz was within [-1 ,1]    
    hxub  = np.where(np.abs(gz)<=1., gx/np.cos(np.arcsin(gz)), 0.)
## bound hx to [-1,1] in case of    
    hx    = np.where(hxub<-1., -1., np.where(hxub>1.,1.,hxub))
    
    gridlatlon[0,:,:] = np.arcsin(gz)*180./np.pi
    gridlatlon[1,:,:] = np.where(gy>0.,np.arccos(hx),-np.arccos(hx))*180./np.pi 
    
    return gridlatlon

block2 = np.zeros(40)
block2[:] = (50, 0, 0, 0, 0, 10, 86, 0, 80, 0, -7020, 0, 0, -7560, 0, 0, 136, 
    7200, 0, 0, 7740, 0, 0, 180, 0, 180, 0, 64, 0, 0, 0, 0, -26500, 0, 0, 
    -44000, 0, 0, 1, 1 )

gridlatlon = get_latlon_shifted(block2)    

