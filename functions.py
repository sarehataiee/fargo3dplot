import numpy as np
import matplotlib.pyplot as plt
from .classes import *

def PlotField(directory, field, nout, fig, ax, cart=False, perturbation=False, residual=False, averaged=False,\
              xlog=False, ylog=False, fieldlog=False, cbar=True, cblabel=True, \
                  settitle=True, **karg):
    """
    This function makes a 2D plot or azimuthally averaged plot of a field. To ease adjucting the 
    plot, the figure and axes must be made before calling this function.

    Parameters
    ----------
    directory : string
        the address where the outputs are located.
    field : string
        name of the field which is going to be plotted. It can be either of 
        "dens" for surface density
        "energy" for energy
        "vx" for azimuthal angular velocity
        "vy" for radial velocity
        "vorticity" for vorticity
        "vortensity" for vorticity divided by density
    nout : int
        the output number
    fig : object
        figure object (it should be made beforehand e.g. by fig, ax = plt.subplots() 
        and is passed to this function)
    ax : object
        the axes this field is going to be demonstrated.
    cart : bool, optional
        If you want to see the plot in xy plane set it to True. The default is False.
    perturbation : bool, optional
        If True, the quantity (field-field_0)/field_0 will be plotted that is the perturbation to the 
        initial output. The default is False.
    residual : bool, optional
        If True, the residual compare to the azimuthally averaged field will be plotted (field-field_av)/field_av.
        The default is False.
    averaged : bool, optional
        If True, the azimuthally averaged (1D) plot will be given. The default is False.
    xlog : bool, optional
        Set to True if you want the x axis to be logarithmic. The default is False.
    ylog : bool, optional
        Set to True if you want the y axis to be logarithmic. The default is False.
    fieldlog : bool, optional
        Set to True if you want the field to be logarithmic. The default is False.
    cbar : bool, optional
        Set to False if you prefer no colorbar in plotting the 2D field
    cblabel : bool, optional
        Set to False if you prefer no colorbar label in plotting the 2D field
    settitle : bool, optional
        Set to False if you prefer no title in plotting the 2D field
    **karg :
        matplotlib plotting arguments either for pcolormesh or plot 

    Returns
    -------
    None.

    """
    grid = Grid(directory)
    r = grid.y
    theta = grid.x
    rmid = grid.ymid
    tmid = grid.xmid
    t2d, r2d = np.meshgrid(theta,r)
    t2dm, r2dm = np.meshgrid(tmid,rmid)
    #
    data = ReadField(directory, field, nout, grid.nx, grid.ny).field
    if fieldlog:
        data = np.log10(data)
        
    cb_label = field
    if cart:
        x = rmid * np.cos(tmid)
        y = rmid * np.sin(tmid)
        X,Y = np.meshgrid(X,Y)
        x_label = "X"
        y_label = "Y"
    else:
        X = tmid
        Y = rmid
        x_label = r"$\theta$"
        y_label = r"$r$"
    if perturbation:
        data_0 = ReadField(directory, field, 0, grid.nx, grid.ny).field
        data = (data-data_0)/data_0
        cb_label += " (perturbation)"
    if residual:
        azi_av_data = np.average(data, axis=1)
        for i in range(grid.ny):
            data[i,:] -= azi_av_data[i]
            data[i,:] /= azi_av_data[i]
        cb_label += " (residual)"
    if averaged:
        azi_av_data = np.average(data, axis=1)
    #
    if not averaged:
        if settitle:
            ax.set_title('{} at output {}'.format(field, nout), fontsize=13)
        cs1 = ax.pcolormesh(X, Y, data, shading='auto', **karg)
        if cbar:
            cb1 = plt.colorbar(cs1, ax=ax, orientation='vertical')
            if cblabel:
                cb1.set_label(cb_label)
    else:
        cs1 = ax.plot(rmid, azi_av_data, label=r"$\rm nout={}$".format(nout), **karg)
        ax.legend()
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    if xlog:
        ax.set_xscale('log')
    if ylog:
        ax.set_yscale('log')        
    
        
        
        
    
