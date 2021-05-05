#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module: units
-------------

This module specifies a container class to provide a set of physical
units (mass, velocity ,length, time) to other modules. The units are
provided by the user and should satisfy the condition that the gravitational
constant is 1 (G = 1.0). 

F. G. Ramon-Fox 2021
Last revision: May 2021
"""

class Units:
    """
    Provides a container for a set of physical units (mass, length, velocity,
    and time) and ensures that they satisfy the condition G=1. This condition
    is checked at instantiation and when a change of units is attempted. The 
    check accepts a small deviation (less than 0.1%) with respect to 1 to allow 
    for truncation errors in user input.
    
    Attributtes
    -----------
    unit_mass : float
        scale unit_mass
    unit_length : float
        scale unit_length
    unit_velocity : float
        scale unit_length
    unit_time : float
        scale unit_time
        
    Usage
    -----
    
    units = Units(umass, ulength, uvel, utime), then
    units.unit_mass returns umass
    
    """
    
    def __init__(self, unit_mass, unit_length, unit_velocity, unit_time):
        
        # Units are validated here during instantiation because after
        # validating type, the units are checked for consistency with G=1.0.
        if not isinstance(unit_mass, float):
            raise TypeError("unit_mass must be a float.")
        if not isinstance(unit_length, float):
            raise TypeError("unit_length must be a float.")
        if not isinstance(unit_velocity, float):
            raise TypeError("unit_velocity must be a float.")
        if not isinstance(unit_time, float):
            raise TypeError("unit_time must be a float.")
        
        if not self.__consistent_units(unit_mass, unit_length, unit_velocity):
            raise ValueError("The entered units are not consistent with G=1.0")
        
        # TODO: Implement a time consistency check.
        
        # We directly assign the attributes here instead of using the
        # property setter since the latter also checks for consistency with G=1.0,
        # so that would cause problems when the object is being instantiated.
        self._unit_mass = unit_mass
        self._unit_length = unit_length
        self._unit_velocity = unit_velocity
        self._unit_time = unit_time
        
    @property
    def unit_mass(self):
        return self._unit_mass
            
    @property
    def unit_length(self):
        return self._unit_length
        
    @property
    def unit_velocity(self):
        return self._unit_velocity
    
    @property
    def unit_time(self):
        return self._unit_time
    
    def update_units(self, new_umass, new_ulength, new_uvel, new_utime):
        """
        Updates the system of units if they are consistent with G=1.0
        
        Parameters
        ----------
        new_umass : float
            new unit mass
        new_ulength : float
            new unit length
        new_uvel : float
            new unit velocity
        new_utime : float
            new unit time
            
        """
        
        if not isinstance(new_umass, float):
            raise TypeError("new_umass must be a float.")
        if not isinstance(new_ulength, float):
            raise TypeError("new_ulength must be a float.")
        if not isinstance(new_uvel, float):
            raise TypeError("new_uvel must be a float.")
        if not isinstance(new_utime, float):
            raise TypeError("new_utime must be a float.")
            
        if not self.__consistent_units(new_umass, new_ulength, new_uvel):
            raise ValueError("The entered units are not consistent with G=1.0")
            
        self._unit_mass = new_umass
        self._unit_length = new_ulength
        self._unit_velocity = new_uvel
        self._unit_time = new_utime
        
    def __consistent_units(self, umass, ulength, uvel):
        """
        Checks if the unit mass, unit length and unit velocity are consistent
        with G=1.0; the check accepts a small deviation (less than 0.1%)
        with respect to 1 to allow for truncation errors in user input.
        
        Parameters
        ----------
        umass : float
            unit mass
        ulength : float
            unit length
        uvel : float
            unit velocity
            
        """
        
        if not isinstance(umass, float):
            raise TypeError("umass must be a float.")
        if not isinstance(ulength, float):
            raise TypeError("ulength must be a float.")
        if not isinstance(uvel, float):
            raise TypeError("uvel must be a float.")
        
        G_ASTRO = 4.301e-6 # kpc/Msun * (km/s)**2
        TOL = 0.001        # Relaxed tolerance to allow for round-off errors.
        
        Gcode = G_ASTRO * (1.0/ulength) * (umass/1.0) * (1.0/uvel)**2
        
        return abs(1.0 - Gcode) < TOL
        
    def __str__(self):
        return str((self.unit_mass, \
                    self.unit_length, \
                    self.unit_velocity, \
                    self.unit_time))