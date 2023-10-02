# -*- coding: utf-8 -*-
"""
Script to combine two colormaps

MADE BY CHRISTIAAN VAN DALUM
"""

def own_colormap(left = 'Blues_r',right = 'Reds', name = 'redblue',totalcolors = 10, paper2 = False, temp = False):
    """ 
    Script to combine two colormaps
    
    MADE BY CHRISTIAAN VAN DALUM
    
    ### INPUT ###
    left: colorset that is on left (DEFAULT = 'Blues_r')
    right: colorset that is on the right (DEFAULT = 'Reds')
    name: name of the new colorset (DEFAULT = 'redblue')
    totalcolors: total number of colors (DEFAULT = 10)
    paper2: Option to use a different setting (DEFAULT = False)
    temp: Option to use a different setting (DEFAULT = False)
    
    ### OUTPUT ###
    newcmp: the new colormap
    
    """

    import numpy as np
    from matplotlib import cm
    from matplotlib.colors import ListedColormap
    
    left_colors = cm.get_cmap(left, 128)
    right_colors = cm.get_cmap(right, 128)
    colors_middle = cm.get_cmap('binary', 128)
    
    newcolors = np.vstack((left_colors(np.linspace(0.0, 1.0, int(totalcolors/2-1))), colors_middle(0.00), colors_middle(0.00),
                           right_colors(np.linspace(0.0, 1.0, int(totalcolors/2-1)))))
    
    if paper2 == True:
        newcolors = np.vstack((left_colors(np.linspace(0.1, 0.9, int(totalcolors/2))),
                               right_colors(np.linspace(0.05, 0.9, int(totalcolors/2)))))
    if temp == True:
        newcolors = np.vstack((left_colors(np.linspace(0.1, 0.85, int(totalcolors/2))), colors_middle(0.00),
                               right_colors(np.linspace(0.15, 0.9, int(totalcolors/2)))))

    newcmp = ListedColormap(newcolors, name=name)
#    test = cm.get_cmap('Greys', 128)
#    newcmp.set_over('k')
#    newcmp.set_under('k')
#    newcmp.set_over(right_colors(1.0))
#    newcmp.set_under(left_colors(0.0))    
    
    test = cm.get_cmap('hot', 128)
    newcmp.set_over(test(0.05))
    test = cm.get_cmap('seismic', 128)
    newcmp.set_under(test(0.0))   
    return newcmp
