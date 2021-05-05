#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module: iofunctions
-------------------

Provides some functions for printing formatted output of the attributes of
several code objects, and some functions for loading data to the code.

F. G. Ramon-Fox 2021
Last Revision: May 2021

"""

from numpy import loadtxt

def print_galaxy_parameters(model):
    """ 
    Display a summary of the parameters of the galaxy model.
    
    Parameters
    ----------
    model : object
        a GalaxyRotationCurve object.
    
    """

    print("  ")
    print("---Halo parameters: ")
    print(" Mh    = " + str(model.Mh))
    print(" c     = " + str(model.c))
    print(" rs    = " + str(model.rs))
    
    print("---Disc parameters: ")
    print(" Md    = " + str(model.Md))
    print(" Rg    = " + str(model.Rg))
    print(" Rd    = " + str(model.Rd))
    print(" gfrac = " + str(model.gfrac))
    print(" sfrac = " + str(model.sfrac))

    print("---Bulge parameters: ")
    print(" Mb    = " + str(model.Mb))
    print(" rb    = " + str(model.rb))
    print("\n")
    
def print_units(units_obj):
    """
    Display the units stored in units_obj.
    
    Parameters
    ----------
    units_obj : object
        a Units object.
        
    """
    
    print("Units: ")
    print(" unit_mass     = " + str(units_obj.unit_mass))
    print(" unit_length   = " + str(units_obj.unit_length))
    print(" unit_velocity = " + str(units_obj.unit_velocity))
    print(" unit_time     = " + str(units_obj.unit_time))
    print("\n")
    
def load_rot_curve_data(path):
    """
    Loads the datapoints of an observed rotation curve stored in path.
    
    the parameter path must include the name of the file appended to the path. 
    In case that the file is in the same directory where the code is being 
    executed, path must include only the file name.
    
    Parameters
    ----------
    path : str
        path, with file name, to the file with the data points
        
    Returns
    -------
    Rdata : ndarray
        radial position (this column should be in kpc in the file).
    vdata : ndarray
        observed rotation curve (this column should be in km/s in the file).
    
    Examples
    --------
    Execution in the same path:
        path="observations.txt"
        Rdata, vdata = load_rot_curve_data(path)
        
    Retrieve data from a different path to the execution directory:
        path="/folder1/folder2/folder3/observations.txt"
        Rdata, vdata = load_rot_curve_data(path)
        
    Notes
    -----
    The file must be structured in the following way:
        # col1: radial pos [kpc] col2: obs vcirc [km/s]
        1.0    50.0
        1.5    100.0
        2.0    200.0
        2.5    220.0
        ...    ...
    
    """
    
    Rdata, vdata = loadtxt(path, unpack=True, usecols=(0,1))
    
    return Rdata, vdata