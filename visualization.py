#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module: visulization
---------------------

Provides several functions for plotting the rotation curve calculated with
the functions provided by a GalaxyRotationCurve object. It also provides
the class PlotData, which serves as a container of plot data and several 
plot parameters that are passed into a function that produces a figure and
writes it to a file.

F. G. Ramon-Fox 2021
Last Revision: May 2021

"""

import numpy as np
import matplotlib.pyplot as plt

class PlotData:
    """
    Class: PlotData
    ---------------
    
    Defines a data container that can be used to store the results of multiple
    radial position vs. circular velocity results for a composite plot. It also
    allows to assign plot configuration parameters such as a label name, marker
    line style, color, X-Y limits, fontsize and legendsize of the plot.
    
    Attributes
    ----------
    n_elements : int
        number of models stored
    Rpos_list : list
        stores the array of radial coordinates of the i-th model.
    vrot_list : list
        stores the array of the circular velocity of the i-th model,
        (vrot_list[i][j] is the circular velocity for Rpos_list[i][j] at
        the j-th position of the individual array).
    label_list : list
        stores the label naming the i-th model.
    marker_list : list
        stores the code or name of the marker to be used for the i-th model.
    ls_list : list
        stores the code or name of the linestyle to be used for the i-th model.
    color_list : list
        stores the name of the color to be used for the i-th model.
    Rpos_lim : list
        pair of limits of the axis of the radial positions (horizontal)
    vrot_lim : list
        pair of limits of the axis of the circular velocity (vertical)
    fontsize : int
        fontsize for axis labels of the plot
    legendsize : int
        fontsize for the legend of the plot
    
    Usage
    -----
    # Instantiate
    pdat = PlotData()
    
    # Assume we have two models
    Rpos = ...
    vgas = ...
    pdat.add_plot_data(Rpos, vgas, label="gas")
    
    Rpos = ...
    vstars = ...
    pdat.add_data(Rpos, vstars, label="stars")
    # See add_plot_data for more details.
    
    plot_composite_rotation_curve(pdata, "figure_name", fmt="png")
    # A file figure_name.png is stored with the figure in the path of execution.
    
    """
    
    def __init__(self):
        
        self._n_elements = 0
        self._Rpos_list = []
        self._vrot_list = []
        self._label_list = []
        self._marker_list = []
        self._ls_list = []
        self._color_list = []
        self._Rpos_lim = []
        self._vrot_lim = []
        self._fontsize = 20
        self._legendsize = 20
        
    @property
    def n_elements(self):
        return self._n_elements
        
    @property
    def Rpos_list(self):
        return self._Rpos_list
    
    @property
    def vrot_list(self):
        return self._vrot_list
    
    @property
    def label_list(self):
        return self._label_list
    
    @property
    def marker_list(self):
        return self._marker_list
    
    @property
    def ls_list(self):
        return self._ls_list
    
    @property
    def color_list(self):
        return self._color_list
    
    @property
    def Rpos_lim(self):
        return self._Rpos_lim
    
    @Rpos_lim.setter
    def Rpos_lim(self, new_Rpos_lim):
        if not isinstance(new_Rpos_lim, list):
            raise TypeError("new_Rpos_lim must be a list.")
        if len(new_Rpos_lim) != 2:
            raise TypeError("new_Rpos_lim must have two elements.")
        
        self._Rpos_lim = new_Rpos_lim
        
    @property
    def vrot_lim(self):
        return self._vrot_lim
    
    @vrot_lim.setter
    def vrot_lim(self, new_vrot_lim):
        if not isinstance(new_vrot_lim, list):
            raise TypeError("new_vrot_lim must be a list.")
        if len(new_vrot_lim) != 2:
            raise TypeError("new_vrot_lim must have two elements.")
            
        self._vrot_lim = new_vrot_lim
    
    @property
    def fontsize(self):
        return self._fontsize
    
    @fontsize.setter
    def fontsize(self, new_fontsize):
        if not isinstance(new_fontsize, int):
            raise TypeError("new_fontsize must be an int.")
        if new_fontsize <= 0:
            raise ValueError("new_fontsize must be a positive int.")
            
        self._fontsize = new_fontsize
    
    @property
    def legendsize(self):
        return self._legendsize
    
    @legendsize.setter
    def legendsize(self, new_legendsize):
        if not isinstance(new_legendsize, int):
            raise TypeError("new_legendsize must be an int.")
        if new_legendsize <= 0:
            raise ValueError("new_legendsize must be a positive int.")
            
        self._legendsize = new_legendsize
    
    def add_plot_data(self, Rpos, vrot, \
                      label=None, marker=None, ls='-', color=None):
        """
        Adds rotation curve data to be plotted and plot configuration parametes
        to the container.
        
        Parameters
        ----------
        Rpos : list or ndarray:
            radial coordinates of the input model.
        vrot : list or ndarray:
            circular velocity of the input model.
            (vrot and Rpos must correspond; vrot is a function of Rpos)
        label : str (optional)
            legend label for this model (default: None)
        marker : str (optional)
            marker label for this model (default: None)
            (use matplotlib markers)
        ls : str (optional)
            key of linestyle to be used for this model (default: '-')
            (use matplotlib styles)
        color : str (optional)
            name of color to be used for this model.
            (use matplotlib color names)
        """
        
        if isinstance(Rpos, list):
            Rpos = np.array(Rpos)
        if not isinstance(Rpos, np.ndarray):
            raise TypeError("Rpos must be a numpy array.")
            
        if isinstance(vrot, list):
            vrot = np.array(vrot)
        if not isinstance(vrot, np.ndarray):
            raise TypeError("vrot must be a numpy array.")
            
        if label is not None and not isinstance(label, str):
            raise TypeError("label must be a string.")
        if marker is not None and not isinstance(marker, str):
            raise TypeError("marker must be a string.")
            
        if not isinstance(ls, str):
            raise TypeError("ls must be a string.")
        if color is not None and not isinstance(color, str):
            raise TypeError("color must be a string.")
            
        
        self._n_elements += 1
        self._Rpos_list.append(Rpos)
        self._vrot_list.append(vrot)
        self._label_list.append(label)
        self._marker_list.append(marker)
        self._ls_list.append(ls)
        self._color_list.append(color)
        
#########

def plot_vrot_vs_radius(Rpos, vrot, fig_name, Rpos_lim=None, vrot_lim=None, \
                        label="", fontsize=20, legendsize=20, fmt="png"):
    """
    Makes a simple plot of circular velocity versus radial position.
    
    Parameters
    ----------
    Rpos : list or ndarray
        radial coordinate array.
    vrot : list or ndarray
        circular velocity array.
    fig_name : str
        base name for the figure output file without extension.
    Rpos_lim : list (optional)
        list with two values: R_lim_min, R_lim_max (default: None)
        if not assigned, the function adapts the limits to the data.
    vrot_lim : None
        list with two values: vrot_lim_min, vrot__lim_max (default: None)
        if not assigned, the function adapts the limits to the data.
    label : str (optional)
        legend label for the plot.    
    fontsize: int (optional)
        defines a fontsize for the axis labels (default: 20 points)
    legendsize: int (optional)
        defines a fontsize for the legend (default: 20 points)
    fmt : str (optional)
        extension of the figure (e.g. png, pdf, eps) (defalut: "png")
        
    Returns
    --------
    filename : str
        a string of the form "fig_name.fmt"
        
    Example
    -------
    set_file_name("figure","jpg") returns "figure.jpg"
    
    """
    
    if isinstance(Rpos, list):
        Rpos = np.array(Rpos, dtype=float)
    if not isinstance(Rpos, np.ndarray):
        raise TypeError("Rpos must be a list or ndarray.")
    if isinstance(vrot, list):
        vrot = np.array(vrot, dtype=float)
    if not isinstance(vrot, np.ndarray):
        raise TypeError("vrot must be a list or ndarray.")
    if not isinstance(fig_name, str):
        raise TypeError("fig_name must be a string.")
    if not isinstance(label, str):
        raise TypeError("label must be a string.")
    if not isinstance(fmt, str):
        raise TypeError("fmt must be a string.")
        
    if label == "":
        label = None
        
    # To use Tex fonts in plot labels and text.
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    
    fig, ax = plt.subplots()
    ax.plot(Rpos, vrot, lw=3, color='black', label=label)
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    
    if Rpos_lim is None:
        Rpos_lim = [0.0, Rpos.max()]
    ax.set_xlim(Rpos_lim)
    
    if vrot_lim is None:
        vrot_lim = [0.0, vrot.max()]    
    ax.set_ylim(vrot_lim)

    ax.set_xlabel(r'$R$ [kpc]', fontsize=fontsize)
    ax.set_ylabel(r'$v_c$ [km/s]', fontsize=fontsize)
    ax.legend(loc=4, fontsize=legendsize)
    
    fig.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
    
    fig_file_name = set_file_name(fig_name, fmt)
    fig.savefig(fig_file_name, format=fmt)    

def plot_composite_rotation_curve(plot_data, fig_name, fmt="png"):
    """
    Returns a file name with the base given by fig_name and the extension
    given by fmt.
    
    Parameters
    ----------
    fig_name : str
        base name string
    fmt : str
        extension of the figure (e.g. png, pdf, eps)
        
    Returns
    --------
    filename : str
        a string of the form "fig_name.fmt"
        
    Example
    -------
    set_file_name("figure","jpg") returns "figure.jpg"
    
    """
    
    if not isinstance(plot_data, PlotData):
        raise TypeError("plot_data must be an instance of PlotData.")
    if not isinstance(fig_name, str):
        raise TypeError("fig_name must be a string.")
    if not isinstance(fmt, str):
        raise TypeError("fmt must be a string.")
    
    # To use TeX fonts in plot labels and text.
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    
    fig, ax = plt.subplots()
    for i in range(plot_data.n_elements):
        ax.plot(plot_data.Rpos_list[i], plot_data.vrot_list[i], \
                label=plot_data.label_list[i], \
                marker=plot_data.marker_list[i], \
                linestyle=plot_data.ls_list[i], 
                lw=3, \
                color=plot_data.color_list[i])
    
    ax.tick_params(axis='both', which='major', labelsize=plot_data.fontsize)
    ax.set_xlim(plot_data.Rpos_lim)
    ax.set_ylim(plot_data.vrot_lim)
    ax.set_xlabel(r'$R$ [kpc]', fontsize=plot_data.fontsize)
    ax.set_ylabel(r'$v_c$ [km/s]', fontsize=plot_data.fontsize)
    ax.legend(loc=4, fontsize=plot_data.legendsize)
    
    fig.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
    
    fig_file_name = set_file_name(fig_name, fmt)
    fig.savefig(fig_file_name, format=fmt)
    
def set_file_name(fig_name, fmt):
    """
    Returns a file name with the base given by fig_name and the extension
    given by fmt.
    
    Parameters
    ----------
    fig_name : str
        base name string
    fmt : str
        extension of the figure (e.g. png, pdf, eps)
        
    Returns
    --------
    filename : str
        a string of the form "fig_name.fmt"
        
    Example
    -------
    set_file_name("figure","jpg") returns "figure.jpg"
    
    """
    
    if not isinstance(fig_name, str):
        raise TypeError("fig_name must be a string.")
    if not isinstance(fmt, str):
        raise TypeError("fmt must be a string.")
    
    if len(fig_name) == 0:
        raise ValueError("fig_name must contain at least one character.")
    if len(fmt) == 0:
        raise ValueError("fmt must contain at least one character.")
    
    filename = fig_name + "." + fmt
    
    return filename
    
    ####################