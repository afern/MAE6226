from matplotlib import pyplot
import numpy as np

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~## PLOTTING FUNCTIONS ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def ezScatter(x, y, x_start, y_start, x_end, y_end):
    """
    x - x values to scatter
    y - y values to scatter
    x_start - x lower limit
    x_end - x upper limit
    y_start - y lower limit
    y_end - y upper limit
    """
    
    width = 10
    height = (y_end - y_start) / (x_end - x_start) * width
    #pyplot.figure(figsize=(width, height))
    pyplot.xlabel('x', fontsize=16)
    pyplot.ylabel('y', fontsize=16)
    pyplot.xlim(x_start, x_end)
    pyplot.ylim(y_start, y_end)
    pyplot.scatter(x, y, color='g', s=80, marker='o')
    return None

def ezPlot(x, y, x_start, y_start, x_end, y_end):
    """
    x - x values to plot
    y - y values to plot
    x_start - x lower limit
    x_end - x upper limit
    y_start - y lower limit
    y_end - y upper limit
    """
    
    width = 10
    height = (y_end - y_start) / (x_end - x_start) * width
    pyplot.figure(figsize=(width, height))
    pyplot.xlabel('x', fontsize=16)
    pyplot.ylabel('y', fontsize=16)
    pyplot.xlim(x_start, x_end)
    pyplot.ylim(y_start, y_end)
    pyplot.plot(x, y)
    return None
    
def ezContourf(x, y, cont, contLevel, x_start, y_start, x_end, y_end):
    """
    x - x values to contour
    y - y values to contour
    cont - contour variable
    contLevel - number of contour lines
    x_start - x lower limit
    x_end - x upper limit
    y_start - y lower limit
    y_end - y upper limit
    """
    width = 10
    height = (y_end - y_start) / (x_end - x_start) * width
    pyplot.figure(figsize=(width, height))
    pyplot.xlabel('x', fontsize=16)
    pyplot.ylabel('y', fontsize=16)
    pyplot.xlim(x_start, x_end)
    pyplot.ylim(y_start, y_end)
    pyplot.contourf(x, y, cont, contLevel)
    return None

def ezStreamline(x, y, U, V, x_start, y_start, x_end, y_end):
    """
    x - x values to scatter
    y - y values to scatter
    U - x velocity
    V - y velocity
    x_start - x lower limit
    x_end - x upper limit
    y_start - y lower limit
    y_end - y upper limit
    """
    width = 10
    height = (y_end - y_start) / (x_end - x_start) * width
    pyplot.figure(figsize=(width, height))
    pyplot.xlabel('x', fontsize=16)
    pyplot.ylabel('y', fontsize=16)
    pyplot.xlim(x_start, x_end)
    pyplot.ylim(y_start, y_end)
    pyplot.streamplot(x, y, U, V, density=2, linewidth=1, arrowsize=1, arrowstyle='->')
    return None

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~## VORTEX FUNCTIONS ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def ezVortexChars(strength, xv, yv, X, Y):
    """
    Returns velocity field generated via vortex, and streamfunction Psi
    *Input params
    -----------
    strength : strength of vortex [FLOAT]
    xv, yv : x, y-coordinate of vortex [FLOAT]
    X, Y : meshed grid of system [2D NP ARRAY]
    
    *Returns
    ---------
    u, v : 2D Numpy arrays of velocity field in x, y-component respectively
    psi : 2D numpy arrays of type [FLOAT] , stream function
    """
    u = strength / (2*np.pi) * (Y - yv) / ((X-xv)**2 + (Y-yv)**2)
    v = -strength / (2*np.pi) * (X - xv) / ((X-xv)**2 + (Y-yv)**2)
    
    psi = strength / (4*np.pi) * np.log((X-xv)**2 + (Y-yv)**2)
    
    return u, v, psi

def complex_potential(strength, a, X, Y):
    """
    Returns the velocity components of a complex potential vortex (infinite sheet)
    
    Parameters:
    ---------
    strength: Strength of the vortex
    a : distance between the vortices
    X, Y: the meshed grid
    
    Returns:
    -------
    u, v : velocity components in cartesian coordinates
    """
    u = +strength/(2*a) * (np.sinh(2*np.pi*Y / a) / (np.cosh(2*np.pi*Y/a) - np.cos(2*np.pi*X/a)))
    v = -strength/(2*a) * (np.sin(2*np.pi*X / a) / (np.cosh(2*np.pi*Y/a) - np.cos(2*np.pi*X/a)))
    
    return u, v