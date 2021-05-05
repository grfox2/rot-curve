#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 
Module: plot_model
------------------- 

Contains the main driver function and some helper functions.

F. G. Ramon-Fox 2021
Last revision: May 2021

"""

import numpy as np
import iofunctions as io
import visualization as vis

from units import Units
from galrotcurve import GalaxyRotationCurve

def main():
    """
    
    This is the main driver function of the code. 
    The user specifies all the relevant parameters in this function 
    (see below.)
    
    """
    
    # Galaxy parameters
    Md = 9.0e9              # disc mass (solar masses)
    Rg = 2.5                # gas disc scale radius (kpc)
    Rd = 2.5                # stellar disc scale radius (kpc)
    gfrac = 0.15            # gas fraction
    sfrac = 0.85            # stellar fraction
    
    # Bulge
    Mb = 3.0e8              # bulge mass (solar masses)
    rb = 0.4                # bulge scale radius (kpc)
    
    # Halo parameters
    Mh = 5.7e11             # halo mass (solar masses)
    c = 4.0                 # concentration
    rs = 33.8               # halo scale radius (kpc)
    
    # Unit parameters.
    unit_mass      = 1.0e5        # Solar Masses
    unit_length    = 0.1          # kpc
    unit_velocity  = 2.074756     # km/s
    unit_time      = 46.96926     # Myr
    
    # Image format
    img_format = "pdf"

    # Include observations data
    include_data = True
    path = "M33_data.txt"
    
    # Radial coordinate parameters
    Rmin = 0.01             # kpc
    Rmax = 15.0             # kpc
    dR = 0.01               # kpc
    
    # Initialize units container, which validates consistency with G = 1.0
    units = Units(unit_mass, unit_length, unit_velocity, unit_time)
    io.print_units(units)

    # Create galaxy model object.
    rcmodel = build_galaxy_model(Md, Rg, Rd, gfrac, sfrac, \
                                 Mb, rb, Mh, c, rs, units)
    io.print_galaxy_parameters(rcmodel)
    
    # plot rotation curve model
    plot_rotation_curve_model(rcmodel, units, Rmin, Rmax, dR=dR, \
                              plot_name="model_curve", fmt=img_format, \
                              include_data=include_data, \
                              data_path=path)

def build_galaxy_model(Md, Rg, Rd, gfrac, sfrac, Mb, rb, Mh, c, rs, units):
    """
    Generates a GalaxyRotationCurve based on the input physical parameters of 
    the model. This parameters are rescaled by the unit system specified in 
    the units object.
    
    Parameters
    ----------
    Md : float or int
        Disc mass in solar masses.
    Rg : float or int
        Gas disc scale radius in kpc.
    Rd : float or int
        Stellar disc scale radius in kpc.
    gfrac : float or int
        Gas fraction (0 to 1.).
    sfrac : float or int
        Stellar fraction (0 to 1.).
    Mb : float or int
        Bulge mass in solar masses.
    rb : float or int
        Bulge scale radius in kpc.
    Mh : float or int
        Dark halo mass in solar masses.
    c : float or int
        Halo concentration parameter.
    rs : float or int
        Halo scale radius in kpc.
    units : object
        Container with the unit system satisfying G=1.
        
    Returns
    -------
    rcmodel : object
        A GalaxyRotationCurve object representing the model.
        
    Usage
    -----
    rcmod = \
        build_galaxy_model(Md, Rg, Rd, gfrac, sfrac, Mb, rb, Mh, c, rs, units)
        
    """
    
    # NOTE: the parameters will be validated at the instantiation
    # GalaxyRotationCurve.
    
    # Create rotationcurve object
    #   Disc & Gas Parameters
    Md = Md/units.unit_mass
    Rg = Rg/units.unit_length
    Rd = Rd/units.unit_length
    
    #   Bulge
    Mb = Mb/units.unit_mass
    rb = rb/units.unit_length
    
    #   Halo parameters
    Mh = Mh/units.unit_mass
    rs = rs/units.unit_length
    
    rcmodel = GalaxyRotationCurve(Md, Rg, Rd, gfrac, sfrac, Mb, rb, Mh, c, rs)
    
    return rcmodel

def plot_rotation_curve_model(rcmod, units, Rmin, Rmax, dR=0.01, \
                              plot_name="rotation_curve", fmt="png", \
                              include_data=False, data_path=None):
    """
    Plots the rotation curve of the model represented by rcmod. It generates
    individual curves of the halo, gas disc, stellar disc, and bulge, as well 
    as the global model. All these results are plotted on the same figure.
    Data points from an observed curve may be optionally included.
    
    Parameters
    ----------
    rcmod : object
        a GalaxyRotationCurve object representing the model.
    units : object
        a Units object, must be the same one used to build rcmod.
    Rmin : float
        minimum radial position to plot.
    Rmax : float
        maximum radial position to plot.
    dR : float (optional)
        separation between radial positions in plot (default: 0.01 in kpc)
    plot_name : str (optional)
        base name of the figure output file, do not include the extension.
        (default: "rotation_curve")
    fmt : str (optional)
        format of the image (e.g. png, pdf, eps)
    include_data : bool (optional)
        if True, reads observed rotation curve data points from data_path.
    data_path : str (optional, necessary if include_data=True)
        filename or path+filename of the observed rotation curve
    
    Example
    -------
    default usage:
        plot_rotation_curve_model(rcmod, units, Rmin, Rmax)
        
    add an observed rotation curve:
        plot_rotation_curve_model(rcmod, units, Rmin, Rmax, '
                include_data=True, data_path="./folder1/folder2/curve.txt")
    
    Notes
    -----
    data_path must point to a two column file with the 1st column containing
    the radial position in kpc, and the second column the rotation curve in km/s.
    See load_rot_curve_data in iofunctions for details.
    
    """
    
    if not isinstance(rcmod, GalaxyRotationCurve):
        raise TypeError("rcmod must be an instance of GalaxyRotationCurve.")
    if not isinstance(units, Units):
        raise TypeError("units must be an instance of Units.")
        
    if not isinstance(Rmin, float) and not isinstance(Rmin, int):
        raise TypeError("Rmin must be a number.")
    if not isinstance(Rmax, float) and not isinstance(Rmax, int):
        raise TypeError("Rmax must be a number.")
    if not isinstance(dR, float):
        raise TypeError("dR must be a float.")
        
    if not isinstance(plot_name, str):
        raise TypeError("plot_name must be a string.")
    if not isinstance(fmt, str):
        raise TypeError("fmt must be a string.")
    if include_data and data_path is None:
        raise TypeError("a data_path must be provided when includ_data is True.")
    if include_data and not isinstance(data_path, str):
        raise TypeError("data_path must be a string.")

    # Generate radial position array        
    R = np.arange(Rmin, Rmax, dR)
    Rcode = R/units.unit_length
    
    # Generate individual curves to visualize the contribution of the 
    # galaxy's components.
    pltdat = vis.PlotData() # Plot data container.
    
    vr_halo = rcmod.get_halo_rotation_curve(Rcode) * units.unit_velocity
    pltdat.add_plot_data(R, vr_halo, label="halo")
    
    vr_gas = rcmod.get_disc_rotation_curve_gas(Rcode) * units.unit_velocity
    pltdat.add_plot_data(R, vr_gas, label="gas disc")
    
    vr_stars = rcmod.get_disc_rotation_curve_stars(Rcode) * units.unit_velocity
    pltdat.add_plot_data(R, vr_stars, label="stellar disc")
    
    vr_bulge = rcmod.get_bulge_rotation_curve(Rcode) * units.unit_velocity
    pltdat.add_plot_data(R, vr_bulge, label="bulge")

    # Get full rotation curve.
    vrot_model = rcmod.get_full_rotation_curve(Rcode) * units.unit_velocity
    pltdat.add_plot_data(R, vrot_model, label="Global", color="black")    
    
    # Load data from observations.
    if include_data:
        Rdata, vdata = io.load_rot_curve_data(data_path)
        pltdat.add_plot_data(Rdata, vdata, \
                             label="observations", ls="none", \
                             marker="o", color="blue")
        
    # Set plot limits and font sizes
    pltdat.Rpos_lim = [0., Rmax]
    pltdat.vrot_lim = [0., 130.]
    pltdat.fontsize = 20
    pltdat.legendeize = 20
        
    # Plot the composite rotation curve.
    vis.plot_composite_rotation_curve(pltdat, "rot_curve_hbd_model", fmt="pdf")
    
    # Plot the simple rotation curve.
    vis.plot_vrot_vs_radius(R, vrot_model, "global_model", label="Global")
    
###########

if __name__ == "__main__":
    main()
