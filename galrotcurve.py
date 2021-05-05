#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" 
Module: galrotcurve
------------------- 

Provides a class to calculate the rotation curve of a galaxy assuming
an axisymmetric model composed by a dark matter halo, disc, and
a spherical bulge.

F. G. Ramon-Fox 2021
Last revision: May 2021

"""

import numpy as np
import scipy.special as sp

class GalaxyRotationCurve:
    """
    Class: GalaxyRotationCurve
    ---------------------------
     
    Provides functions to calculate the rotation curve of a galaxy assuming
    an axisymmetric model composed by a dark matter halo, a stellar disc, 
    a gas disc, and a spherical bulge.
    
    Assumes a natural system of units: G = 1.0, so all the physical quantities
    (position, velocity, mass, etc.) must be scaled by the corresponding
    unit factors. 
    
    For example, Mdisc = Mdisc_solar_masses / unit_mass, where
    Mdisc_solar_masses is the mass of a disc in astrophysical units and Mdisc
    is the mass in natural units, this is the mass to be used in the functions
    of this class. 
     
    Attributes
    ----------
    Md : float or int
        Disc mass in natural units (G=1.0).
    Rg : float or int
        Gas disc scale radius in natural units.
    Rd : float or int
        Stellar disc scale radius in natural units.
    gfrac : float or int
        Gas fraction (0 to 1.)
    sfrac : float or int
        Stellar fraction (0 to 1.)
    Mb : float or int
        Bulge mass in natural units.
    rb : float or int
        Bulge scale radius in natural units.
    Mh : float or int
        Dark halo mass in natural units.
    c : float or int
        Halo concentration parameter.
    rs : float or int
        Halo scale radius in natural units.
    
    Instantiation
    -------------   
     
    mod = GalaxyRotationCurve(Md, Rg, Rd, gfrac, sfrac, Mb, rb, Mh, c, rs)
         
    Notes
    -----
    
    Instantiation and attribute modiciation raise ValueError exceptions when 
    the parameters are out of range, and TypeError when they are not of the 
    correct type.     
    """

    # This class assumes a natural system of units: G = 1.0
    # (it should work for any input parameters in any system 
    # of units satisfying a unitaty G)
    G = 1.0         # Gravitational constant
    
    def __init__(self, Md, Rg, Rd, gfrac, sfrac, Mb, rb, Mh, c, rs):
        
        # Disc parameters
        self.Md = Md
        self.Rg = Rg
        self.Rd = Rd
        self.gfrac = gfrac
        self.sfrac = sfrac
        
        # Bulge parameters
        self.Mb = Mb
        self.rb = rb
        
        # Halo parameters
        self.Mh = Mh
        self.rs = rs    
        self.c = c
        
    # Disc property getters/setters
    @property
    def Md(self):
        return self._Md
    
    @Md.setter
    def Md(self, new_Md):
        if new_Md < 0.:
            raise ValueError("new_Md must be positive.")
        if not isinstance(new_Md, float) and not isinstance(new_Md, int):
            raise TypeError("new_Md must be a number.")    
            
        self._Md = float(new_Md)
        
    @property
    def Rg(self):
        return self._Rg
    
    @Rg.setter
    def Rg(self, new_Rg):
        if new_Rg < 0.:
            raise ValueError("new_Rg must be positive.")
        if not isinstance(new_Rg, float) and not isinstance(new_Rg, int):
            raise TypeError("new_Rg must be a number.")
            
        self._Rg = float(new_Rg)
        
    @property
    def Rd(self):
        return self._Rd
    
    @Rd.setter
    def Rd(self, new_Rd):
        if new_Rd < 0.:
            raise ValueError("new_Rd must be positive.")
        if not isinstance(new_Rd, float) and not isinstance(new_Rd, int):
            raise TypeError("new_Rd must be a number.")
        
        self._Rd = float(new_Rd)
        
    @property
    def gfrac(self):
        return self._gfrac
    
    @gfrac.setter
    def gfrac(self, new_gfrac):
        if new_gfrac < 0. or new_gfrac > 1.:
            raise ValueError("new_gfrac must be between 0 and 1 inclusive.")
        if not isinstance(new_gfrac, float):
            raise TypeError("new_gfrac must be a float.")
        
        self._gfrac = float(new_gfrac)
        
    @property
    def sfrac(self):
        return self._sfrac
    
    @sfrac.setter
    def sfrac(self, new_sfrac):
        if new_sfrac < 0. or new_sfrac > 1.:
           raise ValueError("new_sfrac must be between 0 and 1 inclusive.") 
        if not isinstance(new_sfrac, float):
            raise TypeError("new_sfrac must be a float.")
        
        self._sfrac = float(new_sfrac)
        
    # Bulge properties getters/setters
    @property
    def Mb(self):
        return self._Mb
    
    @Mb.setter
    def Mb(self, new_Mb):
        if new_Mb < 0.:
            raise ValueError("new_Mb must be positive.")  
        if not isinstance(new_Mb, float) and not isinstance(new_Mb, int):
            raise TypeError("new_Mb must be a number.")
            
        self._Mb = float(new_Mb)
        
    @property
    def rb(self):
        return self._rb
    
    @rb.setter
    def rb(self, new_rb):
        if new_rb < 0.:
            raise ValueError("new_rb must be positive.")
        if not isinstance(new_rb, float) and not isinstance(new_rb, int):
            raise TypeError("new_rb must be a number.")
            
        self._rb = float(new_rb)
        
    # Halo peroperties getters/setters.
    @property
    def Mh(self):
        return self._Mh
    
    @Mh.setter
    def Mh(self, new_Mh):
        if new_Mh < 0.:
            raise ValueError("new_Mh must be positive.")
        if not isinstance(new_Mh, float) and not isinstance(new_Mh, int):
            raise TypeError("new_Mh must be a number.")
            
        self._Mh = float(new_Mh)
        
    @property
    def rs(self):
        return self._rs
    
    @rs.setter
    def rs(self, new_rs):
        if new_rs < 0.:
            raise ValueError("new_rs must be positive.")
        if not isinstance(new_rs, float) and not isinstance(new_rs, int):
            raise TypeError("new_rs must be a number.")
            
        self._rs = float(new_rs)
        
    @property
    def c(self):
        return self._c
    
    @c.setter
    def c(self, new_c):
        if new_c < 0.:
            raise ValueError("new_c must be positive.") 
        if not isinstance(new_c, float) and not isinstance(new_c, int):
            raise TypeError("new_c must be a number.")
            
        self._c = float(new_c)
        
    def get_halo_rotation_curve(self, r):
        """ 
        Returns the rotation curve of the dark matter halo as a function
        of r, assuming a Navarro-Frenk-White (1996) profile and spherical 
        symmetry.
        
        Parameters
        ----------
        r : array-like
            radial position array in natural (code) units.
            
        Returns
        -------
        vh_rot :
            rotation curve of the dark halo in natural units.
            
        Notes
        -----
        Multiply vh_rot by your unit velocity to convert to physical units.
        """
        
        if not isinstance(r, np.ndarray):
            raise TypeError("r must be a ndarray.")
        
        x = r/self._rs
        gc = np.log(1.0 + self._c) - self._c/(1.0 + self._c)
        vh_rotsqd = (GalaxyRotationCurve.G * self._Mh/self._rs) * 1.0/gc \
            * ( np.log(1.0+x) - x/(1.0+x) )/x
        vh_rot = np.sqrt(vh_rotsqd)
        
        return vh_rot
        
    def get_disc_rotation_curve_gas(self, R):
        """ 
        Returns the rotation curve of a gas disc with an exponential surface
        density profile. The mass of the gas disc is Mgas = Md * gasfrac.
        
        Parameters
        ----------
        R : array-like
            radial position array in natural (code) units.
            
        Returns
        -------
        vd_rot :
            rotation curve of the gas disc in natural units.
            
        Notes
        -----
        Multiply vd_rot by your unit velocity to convert to physical units.
        """

        if not isinstance(R, np.ndarray):
            raise TypeError("R must be a ndarray.")

        Mgas = self._Md * self._gfrac
        vd_rotsqd = self.__get_disc_rotation_curve(Mgas, self._Rd, R)
        vd_rot = np.sqrt(vd_rotsqd)
        
        return vd_rot
        
    def get_disc_rotation_curve_stars(self, R):
        """ 
        Returns the rotation curve of a stellar disc with an exponential surface
        density profile. The mass of the stellar disc is Mstars = Md * gasfrac.
        
        Parameters
        ----------
        R : array-like
            radial position array in natural (code) units.
            
        Returns
        -------
        vd_rot :
            rotation curve of the stellar disc in natural units.
            
        Notes
        -----
        Multiply vd_rot by your unit velocity to convert to physical units.
        """

        if not isinstance(R, np.ndarray):
            raise TypeError("R must be a ndarray.")

        Mstars = self._Md * self._sfrac
        vd_rotsqd = self.__get_disc_rotation_curve(Mstars, self._Rd, R)        
        vd_rot = np.sqrt(vd_rotsqd)
        
        return vd_rot
    
    def __get_disc_rotation_curve(self, Mdisc, Rdisc, R):
        """
        Returns the squared of the circular velocity of a disc with an 
        exponential surface density profile.
        
        Parameters
        ----------
        R : array-like
            radial position array in natural (code) units.
            
        Returns
        -------
        vd_rotsqd :
            squared of the circular velocity of the disc in natural units.
            
        Notes
        -----
        Do not forget to apply the square root and multiply by velocity unit.
        This function is aimed for internal use only.
        """
        
        if not isinstance(Mdisc, float):
            raise TypeError("Mdisc must be a float.")
        if not isinstance(Rdisc, float):
            raise TypeError("Rdisc must be a float.")
        if not isinstance(R, np.ndarray):
            raise TypeError("R must be a ndarray.")
        
        y = R/(2.0 * Rdisc)
        vd_rotsqd = (2.0 * GalaxyRotationCurve.G * Mdisc/Rdisc) * y**2.0 \
            * (sp.iv(0,y)*sp.kv(0,y) - sp.iv(1,y)*sp.kv(1,y))
        
        return vd_rotsqd
    
    def get_bulge_rotation_curve(self, r):
        """ 
        Returns the rotation curve of a spherical bulge with the density
        profile of Hernquist (1990).
        
        Parameters
        ----------
        r : array-like
            radial position array in natural (code) units.
            
        Returns
        -------
        vb_rot :
            rotation curve of the bulge in natural units.
            
        Notes
        -----
        Multiply vb_rot by your unit velocity to convert to physical units.
        """
        
        if not isinstance(r, np.ndarray):
            raise TypeError("r must be a ndarray.")
        
        xi = r/self._rb
        vb_rotsqd = (GalaxyRotationCurve.G * self._Mb/self._rb) * \
            (xi/(1.0 + xi)**2.0)
        vb_rot = np.sqrt(vb_rotsqd)
        
        return vb_rot

    def get_full_rotation_curve(self, r):
        """ 
        Returns the total rotation curve of the specified galactic model.
        
        Parameters
        ----------
        r : array-like
            radial position array in natural (code) units.
            
        Returns
        -------
        vrot :
            rotation curve of the galactic model in natural units..
            
        Notes
        -----
        Multiply vrot by your unit velocity to convert to physical units.
        """
        
        if not isinstance(r, np.ndarray):
            raise TypeError("r must be a ndarray.")
        
        vh_rotsqd = (self.get_halo_rotation_curve(r))**2
        vd_rotsqd = (self.get_disc_rotation_curve_stars(r))**2.0
        vg_rotsqd = (self.get_disc_rotation_curve_gas(r))**2.0
        vb_rotsqd = (self.get_bulge_rotation_curve(r))**2.0
        
        return np.sqrt(vh_rotsqd + vd_rotsqd + vg_rotsqd + vb_rotsqd)
    
    def __str__(self):
        return "({0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9})".\
            format(self._Md, self._Rg, self._Rd, self._gfrac, self._sfrac, \
                   self._Mb, self._rb, self._Mh, self._c, self._rs)
        
# End GalaxyRotationCurve
########